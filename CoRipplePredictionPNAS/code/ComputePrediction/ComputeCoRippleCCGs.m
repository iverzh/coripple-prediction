
%% ComputeCoRippleCCGs.m
% compute cross correlograms between all units during co-ripple and control periods

clear
close all
clc

addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))

% parpool(12)
%%
subject = 'E1';
state = 'NREM'; %NREM or wake
aList = 1:nUnits;

outputDir = sprintf('./out/%s/CCG/%s', subject, state);
if ~isfolder(outputDir); mkdir(outputDir); end

units = load(); %see readMe for instructions on how to format unit matrix

matExportFolder = ''; %folder for output of ripple detection.
load(fullfile(matExportFolder,filename));

load('') %load stage mask

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

rippleChannels = cellfun(@(x) str2double(x), rippleStats.chanLabels);

%%
tic 
computeCoRipPredict = true; %compute co-ripple and no-ripple CCGs
computeSingleRipPredict = true; %compute CCG when only A unit has a ripple
cd(outputDir)
close all
win = [-750 0];
ripp = [];
rippA = [];
null   = [];

fname = sprintf('%s_%i-%i_%s.log',subject,aList(1),aList(end), state);
fileID = fopen(fullfile(outputDir,fname), 'w');
fprintf('Computing Spike CCGs\n');
fprintf(num2str(aList))
for ai = 1:length(aList)
    a = aList(ai);
    ripp = [];
    rippA = [];
    null   = [];    %% 
    ua = unitsAll(a);
    
    chA = units{ua,1};
    map = chMap == chA;
    [Xaii, Yaii] = find(map);
    rMaska = rippMask(rippleChannels == chA,:);

    for b = 1:length(unitsAll)

        if a ~= b
            clc

            t = datetime('now');
            s = datestr(t);
            s(s == ' ') = '-';


            fprintf(fileID, '%s \n', s);
            fprintf(fileID, 'chan A %i out of %i \n', a, length(unitsAll));
            fprintf(fileID, 'chan B %i out of %i \n\n', b, length(unitsAll));


            ub = unitsAll(b);
            chB = units{ub,1};
            
            map = chMap == chB;
            rMaskb = logical(rippMask(rippleChannels == chB,:));
            if computeCoRipPredict
                rMaskab = rMaska & rMaskb;
    
                APa = units{ua,2} * rippleStats.fs;
                APb = units{ub,2} * rippleStats.fs;
                APa(APa > rippleStats.recordingLength) = [];
                APb(APb > rippleStats.recordingLength) = [];
    
                APbr = APb(StageMask_1kHz(round(APb)));
                keepNull = ~rMaskb(round(APb)) & StageMask_1kHz(round(APb));
                APbn = APb(keepNull);
    
                APar = APa(logical(rMaskab(round(APa)))); %find spikes in overlapping ripples. 
    
                APan_ii = find(~rMaska(round(APa)) & StageMask_1kHz(round(APa))); %find spikes in chA . 
     
                if length(APan_ii)/length(APar) > 30 
                    iiKeep = randsample(length(APan_ii)-1,30*length(APar));
                else
                    iiKeep = 1:length(APan_ii);
                end
                
                APan = APa(APan_ii(sort(iiKeep)));
    
                spikeHistoryRipple  = nan(1e6,length(win(1):win(2))+1);
                cRipple = 1;
    
                spikeHistoryControl = nan(1e6,length(win(1):win(2))+1);
                cControl = 1;
                
                spikeRBbRipple  = zeros(1e6,length(win(1):win(2)));
                spikeRBaRipple  = zeros(1e6,length(win(1):win(2)));
                spikeRBbControl = zeros(1e6,length(win(1):win(2)));
                spikeRBaControl = zeros(1e6,length(win(1):win(2)));
    
                for ii = 1:length(APar)
                  
                   checkStage = StageMask_1kHz(round(APar(ii)));
               
    
                   if checkStage 
                       checkAPb = APbr > (APar(ii) + win(1)) & APbr < (APar(ii) + win(2));
                       if any(rMaskb(round(APbr(checkAPb))))
                           checkAPb = APbr(checkAPb) - APar(ii);
                           checkAPb(checkAPb == 0) = [];
        
        
        
                           spikeHistory = zeros(size(win(1):win(2)));
                           spikeHistory(round(checkAPb-win(1)) + 1) = 1;
                           if APar(ii)+win(1) > 0
                               spikeHistoryRipple(cRipple,:) = [APar(ii) spikeHistory];
       
                           else
                               spikeHistoryRipple(cRipple,:) = [nan nan(size(win(1):win(2)))];

                           end
                           cRipple = cRipple + 1;
                       end
                   end
                end
    
                for ii = 1:length(APan)
                   checkStage = StageMask_1kHz(round(APan(ii)));
    
    
                   if checkStage
                       checkAPb = APbn > (APan(ii) + win(1)) & APbn < (APan(ii) + win(2));
                       checkAPb = APbn(checkAPb) - APan(ii);
                       checkAPb(checkAPb == 0) = [];               
    
                       spikeHistory = zeros(size(win(1):win(2)));
                       spikeHistory(round(checkAPb-win(1)) + 1) = 1;
                       if APan(ii)+win(1) > 0
                           spikeHistoryControl(cControl,:) = [APan(ii) spikeHistory];
    
                       else
                           spikeHistoryControl(cControl,:) = [nan nan(size(win(1):win(2)))];
                  
                       end
    
                       cControl = cControl + 1;
    
                   end            
                end

                spikeHistoryRipple(cRipple:end,:) = [];
                spikeRBaRipple(cRipple:end,:) = [];
                spikeRBbRipple(cRipple:end,:) = [];
                spikeHistoryControl(cControl:end,:) = [];
                spikeRBaControl(cControl:end,:) = [];
                spikeRBbControl(cControl:end,:) = [];
    
    
    
                ripp.History{a,b} = spikeHistoryRipple;
                null.History{a,b} = spikeHistoryControl;
           

                
            end

            if computeSingleRipPredict
                rMaskaNob = rMaska & ~rMaskb;

    
                APa = units{ua,2} * rippleStats.fs;
                APb = units{ub,2} * rippleStats.fs;
                APa(APa > rippleStats.recordingLength) = [];
                APb(APb > rippleStats.recordingLength) = [];
    
                keepNull = ~rMaskb(round(APb)) & StageMask_1kHz(round(APb));
                APbn = APb(keepNull);
    
                APar = APa(logical(rMaskaNob(round(APa)))); %find spikes with ripple in a and no ripple in b. 
           
                spikeHistoryRipple  = nan(1e6,length(win(1):win(2))+1);
                cRipple = 1;

                for ii = 1:length(APar)
                  
                   checkStage = StageMask_1kHz(round(APar(ii)));
               
    
                   if checkStage 
                       checkAPb = APbn > (APar(ii) + win(1)) & APbn < (APar(ii) + win(2));
                       if all(~rMaskb(round(APbn(checkAPb)))) && ~isempty(checkAPb)
                           checkAPb = APbn(checkAPb) - APar(ii);
                           checkAPb(checkAPb == 0) = [];
        
        
        
                           spikeHistory = zeros(size(win(1):win(2)));
                           spikeHistory(round(checkAPb-win(1)) + 1) = 1;
                           if APar(ii)+win(1) > 0

                               spikeHistoryRipple(cRipple,:) = [APar(ii) spikeHistory];
        
                           else
                               spikeHistoryRipple(cRipple,:) = [nan nan(size(win(1):win(2)))];

                           end
                           cRipple = cRipple + 1;
                       end
                   end
                   
                end

                spikeHistoryRipple(cRipple:end,:) = [];
                rippA.History{a,b} = spikeHistoryRipple;

            end


        end
        
          
    end
    
    if computeCoRipPredict
        nRipp = [];
        nNull = [];
        for bb = 1:size(ripp.History,2)
            if bb ~= a
                nRipp = [nRipp, size(ripp.History{a, bb},1)];
                nNull = [nNull, size(null.History{a, bb},1)];
            end
    
        end



        fname = sprintf('%s_CCGs_%s_unit_%i.mat',subject, state, a);
        save(fullfile(outputDir,fname), 'null', 'ripp', 'nRipp','nNull','-v7.3')
        clear null; clear ripp;
    end
    if computeSingleRipPredict
        nRippA = [];
        for bb = 1:size(rippA.History,2)
            if bb ~= a
                nRippA = [nRippA, size(rippA.History{a, bb},1)];
            end
    
        end

        fname = sprintf('%s_CCGs_%s_unit_%i.mat',subject, state, a);
        save(fullfile(outputDir,fname), 'rippA', 'nRippA','-append')
        clear rippA;
    end
    
end

toc