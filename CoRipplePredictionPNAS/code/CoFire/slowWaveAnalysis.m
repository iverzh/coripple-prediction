


close all
clear
clc
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/RippleDetection'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))

%% Inputs

subject = 'T11'; %'T11';
state = 'NREM';
rawDataPath = sprintf('/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/%s/data_1kHz/', subject);
filename = sprintf('%s_1kHz_unitsRemoved_NREM.mat', subject);

exportSO = sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/SlowWaves/%s', subject);
if ~isfolder(exportSO); mkdir(exportSO); end

runSODetect = false;
brainGate = true; %true;
arrayConfig = 'single'; %'dual';

if ~brainGate
    units = LoadSpikeTimes(subject,'CellExplorer','medial');
    
    matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
    load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMOnly.mat']));
    rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
    for chRipp = 1:size(rippMask,1) 
        if ~isempty(rippleStats.window{chRipp})
            iS = rippleStats.window{chRipp}(:,1);
            iE = rippleStats.window{chRipp}(:,2);
        
            for ii = 1:length(iE)
                rippMask(chRipp,iS(ii):iE(ii)) = 1;
            end
        end
        
    end
    
    
    if strcmp(subject,'MG63')
           StageMask_1kHz(4229990:15900000) = 0;
    end
else
    if strcmp(arrayConfig ,'dual')
        unitsM = LoadSpikeTimes(subject,'CellExplorer','medial');
        unitsL = LoadSpikeTimes(subject,'CellExplorer','lateral');
        
        for ii = 1:size(unitsL,1); unitsL{ii,1} = unitsL{ii,1} + 96; end
        
        units = [unitsM; unitsL];
    
        matExportFolder = '/space/seh9/2/halgdev/projects/iverzh/ripples/matFiles';
        rippMedial   = load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMonly_medial.mat']));
        rippLateral  = load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMonly_lateral.mat']));
        
        rippMaskM = zeros(length(rippMedial.rippleStats.chanLabels), rippMedial.rippleStats.recordingLength);
        rippMaskL = zeros(length(rippLateral.rippleStats.chanLabels), rippLateral.rippleStats.recordingLength);
        for chRipp = 1:size(rippMaskM,1) 
            if ~isempty(rippMedial.rippleStats.window{chRipp})
                iS = rippMedial.rippleStats.window{chRipp}(:,1);
                iE = rippMedial.rippleStats.window{chRipp}(:,2);
            
                for ii = 1:length(iE)
                    rippMaskM(chRipp,iS(ii):iE(ii)) = 1;
                end
            end
        
            if ~isempty(rippLateral.rippleStats.window{chRipp})
                iS = rippLateral.rippleStats.window{chRipp}(:,1);
                iE = rippLateral.rippleStats.window{chRipp}(:,2);
            
                for ii = 1:length(iE)
                    rippMaskL(chRipp,iS(ii):iE(ii)) = 1;
                end
            end
            
        end
        
        rippMask = [rippMaskM; rippMaskL];
    elseif strcmp(arrayConfig, 'single')
        units = LoadSpikeTimes(subject,'CellExplorer','medial');
        matExportFolder = '/space/seh9/2/halgdev/projects/iverzh/ripples/matFiles';
        load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMonly_medial.mat']));
        rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
        for chRipp = 1:size(rippMask,1) 
            if ~isempty(rippleStats.window{chRipp})
                iS = rippleStats.window{chRipp}(:,1);
                iE = rippleStats.window{chRipp}(:,2);
            
                for ii = 1:length(iE)
                    rippMask(chRipp,iS(ii):iE(ii)) = 1;
                end
            end
            
        end

    end

   
    
end
load(sprintf('/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/%s/data_1kHz/%s_stageMask_%s.mat', subject, subject, state))

rippMask(:,~StageMask_1kHz) = 0;
%% Run SO wave detection
USwin = 150; 
if runSODetect
    load(fullfile(rawDataPath,filename))
    data(isnan(data)) = 0;
    [ SOmat, summary ] = detectSO(-data, 1e3, [], 20);
    save(fullfile(exportSO, 'SO_NREM.mat'), "SOmat", "summary", '-v7.3')
end

load(fullfile(exportSO, 'SO_NREM.mat'));

%%
USmask = zeros(size(rippMask,1), size(rippMask,2));
for ch = 1:size(USmask,1) 
    SOch = SOmat(SOmat(:,3) == ch, :);
    USch = SOch(SOch(:,2) > 0, 1);
    if ~isempty(USch)
        for ii = 1:length(USch)
            USmask(ch,USch(ii)-USwin:USch(ii)) = 1;
        end
    end
    
end
USmask(:,~StageMask_1kHz) = 0;

%%
singleUnits = cellfun(@(X) any(strcmp(X, {'pyr', 'int'})), units(:,3))
Usall = sum(USmask);
b = mask2bounds(Usall > 32);

for iB = 1:length(b)
    

end




%%
USrippCoFire = nan(length(units), length(units));
ripprippCoFire = nan(length(units), length(units));
USUSCoFire = nan(length(units), length(units));
nullCoFire = nan(length(units), length(units));
coFireThresh = 10;

reverseStr = '';
for iiA = 1:length(units)
    typeA = units{iiA,3};
    chA = units{iiA,1};
    USa = USmask(chA,:);
    rippA = rippMask(chA,:);
    timesA = round(units{iiA,2} * 1e3); %multiply by sampling rate
    if strcmp(typeA,'pyr') || strcmp(typeA,'int') && chA <= 96
        for iiB = 1:length(units)
            typeB = units{iiB,3};
            chB = units{iiB,1};
            USb = USmask(chB,:);
            rippB = rippMask(chB,:);
            timesB = round(units{iiB,2} * 1e3); %multiply by sampling rate

            if (strcmp(typeB,'pyr') || strcmp(typeB,'int')) && chB <= 96
                USripp = USb & USa & rippA & rippB;
                rippripp = ~USb & ~USa & rippA & rippB;
                USUS = USb & USa & ~rippA & ~rippB;
                null = ~USb & ~USa & ~rippA & ~rippB & StageMask_1kHz;
                
                cond = USripp;
                eventA = timesA(cond(timesA));
                eventB = timesB(cond(timesB));
                s = eventA - eventB';
                coFire = sum(s(:) <= coFireThresh & s(:) >= 0);
                USrippCoFire(iiA,iiB) = coFire / sum(cond);

                cond = rippripp;
                eventA = timesA(cond(timesA));
                eventB = timesB(cond(timesB));
                s = eventA - eventB';
                coFire = sum(s(:) <= coFireThresh & s(:) >= 0);
                ripprippCoFire(iiA,iiB) = coFire / sum(cond);

                cond = USUS;
                eventA = timesA(cond(timesA));
                eventB = timesB(cond(timesB));
                s = eventA - eventB';
                coFire = sum(s(:) <= coFireThresh & s(:) >= 0);
                USUSCoFire(iiA,iiB) = coFire / sum(cond);

                cond = null;
                eventA = timesA(cond(timesA)); 
                eventB = timesB(cond(timesB));
                try
                    s = eventA - eventB';
                    coFire = sum(s(:) <= coFireThresh & s(:) >= 0);
                catch
                    warning('co fire matrix is too large. Using first night only')
                    eventA(eventA > 14413800) = [];
                    eventB(eventB > 14413800) = [];
                    s = eventA - eventB';
                    coFire = sum(s(:) <= coFireThresh & s(:) >= 0);
                    cond(14413801:end) = 0;
                end
                nullCoFire(iiA,iiB) = coFire / sum(cond);

               

            end


        end

    

    end
    percentDone = 100 * iiA / length(units);
    msg = sprintf('Running US ripple co fire analysis: %3.1f \n', percentDone); %Don't forget this semicolon
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));

end
filename = sprintf('%s_SO_CoFire.mat', [subject, '-', arrayConfig]);
save(fullfile(exportSO,filename), 'nullCoFire', "USUSCoFire", "ripprippCoFire", "USrippCoFire", '-v7.3')

%%

spikeMaskPyr = zeros(size(rippMask));
spikeMaskInt = zeros(size(rippMask));
spikeMask = zeros(size(rippMask));

for uii = 1:length(units)
    uType = units{uii,3};
    ch = units{uii,1};
    spikes = round(units{uii,2} * 1e3);
%     spikes(spikes < 15900000) = [];
    if strcmp(uType,'pyr')
        spikeMaskPyr(ch,spikes) = 1;
        secs = (rippleStats.recordingLength - 15900000) / 1e3;
        FR = length(spikes) / secs;
        sprintf('%s %.4f Hz', uType,FR)
        spikeMask(ch,spikes) = 1;

    elseif strcmp(uType,'int')

        spikeMaskInt(ch,spikes) = 1;

        secs = (rippleStats.recordingLength - 15900000) / 1e3;
        FR = length(spikes) / secs;
        sprintf('%s %.4f Hz', uType,FR)
        spikeMask(ch,spikes) = 1;

    end

end



US_noRipp = logical(USmask) & logical(rippMask);
% US_noRipp =  logical(rippMask);
null = ~USmask & ~rippMask & StageMask_1kHz';
%%

for ch = 1:96
    b = mask2bounds(US_noRipp(ch,:));
    dur = b(:,2) - b(:,1);
%     b(dur < 150,:) = [];
    US_noRipp(ch,:) = bounds2mask(b, length(US_noRipp));
    
    FR = sum(spikeMask(ch,US_noRipp(ch,:))) / sum(US_noRipp(ch,:)) * 1e3;
    fprintf('%.2f  ', FR)
    FR = sum(spikeMask(ch,null(ch,:))) / sum(null(ch,:)) * 1e3;
    fprintf('%.4f  \n', FR)

end


%%
figure;
for iB = 1:length(b)
%     loc = SOch(iB,1);
%     mag = SOch(iB,2);
%     if mag > 0 && any(rippMask(ch, loc-150:loc))
    plot(-data(ch, b(iB,1)-200:b(iB,2)+200)); hold on;
%         plot(-data(ch, loc-200:loc+200)); hold on;
        vline(350)
    

        waitforbuttonpress; clf;

%     end
end







