clear
clc
close all


% parpool(12)
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/ripple-detection/code/utils'))
%%
subject = 'T11';
state = 'wake';
units = LoadSpikeTimes(subject,'CellExplorer','medial');
matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
% load(fullfile(matExportFolder,[subject,'_ripple_stats_wake_NC_wakeOnly.mat']));
% load(fullfile(matExportFolder,[subject,'_ripple_stats_wake_NC_wakeOnly.mat    ']));
% load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMonly_medial.mat']));
% load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMonly_lateral.mat']));
load(fullfile(matExportFolder,[subject,'_ripple_stats_wake_NC_wakeonly_medial.mat']));
% load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMOnly.mat']));
% load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMOnly_HGevents.mat']));

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
%         secs = (rippleStats.recordingLength - 15900000) / 1e3;
%         FR = length(spikes) / secs;
        sprintf('%s %.4f Hz', uType)
        spikeMask(ch,spikes) = 1;

    elseif strcmp(uType,'int')

        spikeMaskInt(ch,spikes) = 1;

%         secs = (rippleStats.recordingLength - 15900000) / 1e3;
%         FR = length(spikes) / secs;
        sprintf('%s %.4f Hz', uType)
        spikeMask(ch,spikes) = 1;

    end

end

load(sprintf('/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/%s/data_1kHz/%s_stageMask_%s.mat', subject, subject, state))


if strcmp(subject,'MG63')
       StageMask_1kHz(4229990:15900000) = 0;
end

rippMask(:,~StageMask_1kHz) = 0;
%%
% coFireThresh = [1, 5, 11, 15, 22];
coFireProbCoR = cell(length(units),length(units));
coFireProbCoR_Ana = cell(length(units),length(units));
coFireProbCoR_Shuff = cell(length(units),length(units));
coRippleFRa = nan(length(units),length(units));
coRippleFRb = nan(length(units),length(units));
coRippleDur = nan(length(units),length(units));
% coFireProbNoR = nan(length(units),length(units));
% coFireProbNoR_Ana = nan(length(units),length(units));
interactionType = cell(length(units),length(units));
c = 0;
for uAi = 1:length(units)
    for uBi = 1:length(units)
        clc
        c = c+1;
        fprintf('co fire analysis \n ==== %.4f ====\n' ,c/length(units)^2)
%         uA = unitsA(uAi); uB = unitsA(uBi);
        uA = uAi; uB = uBi;
        chA = units{uA,1}; chB = units{uB,1};
        typeA = units{uA,3}; typeB = units{uB,3};
        
        if chA ~= chB && any(strcmp(typeA,{'pyr', 'int'})) && any(strcmp(typeB,{'pyr', 'int'}))
%         if chA ~= chB && any(strcmp(typeA,{'pyr'})) && any(strcmp(typeB,{'pyr'}))
            uTa = units{uA,2} * 1e3;
            uTb = units{uB,2} * 1e3;

            coRipMaskO = rippMask(chA, :) & rippMask(chB, :);
            b = mask2bounds(coRipMaskO);
            coRipMask = bounds2mask(b,length(coRipMaskO),100000)';

            uTcoRa = uTa(coRipMask(round(uTa)));
            uTcoRb = uTb(coRipMask(round(uTb)));

            if ~isempty(uTcoRa) && ~isempty(uTcoRb)
                
                
                coRipFRa = sum(coRipMask(round(uTa))) / sum(coRipMask) * 1e3; %Hz
                coRipFRb = sum(coRipMask(round(uTb))) / sum(coRipMask) * 1e3; %Hz

                s = uTcoRa - uTcoRb';
%                 s = abs(s);
%                 for iThresh = 1:length(coFireThresh)
%                 coF = sum(s(:) <= coFireThresh & s(:) >= 0);

%                 coFire = coF / sum(coRipMask) * 1e3;
                coFireChanceAna = 2 * (coRipFRa * coRipFRb) * [1:50] /1e3;
                coFireChance = computeChanceCoFire(uTcoRa, uTcoRb, coRipMask, [], s, 1000);
                
                N = histcounts(s, -50.5:50.5);
                coFireProbCoR{uA,uB} = N;
                coFireProbCoR_Shuff{uA,uB} = coFireChance; %(mean(coFireChance) / sum(coRipMask) * 1e3);
                coFireProbCoR_Ana{uA,uB} = coFireChanceAna;
                coRippleFRa(uA,uB) = coRipFRa;
                coRippleFRb(uA,uB) = coRipFRb;
                coRippleDur(uA,uB) = sum(coRipMask);

                coShuff = sum(coFireChance(:,51-coF:51+coF), 2);
                figure;
                histogram(coShuff); hold on;
                vline(sum(N(51-coF:51+coF))); hold on;
                waitforbuttonpress; close;  

%                 if strcmp(typeA, 'pyr') && strcmp(typeB, 'int') && coF > 10
%                      
%                     uA
%                     uB
%                     coF
%                     pause
%                 end
    
%                 noRipMask = ~rippMask(chA, :) & ~rippMask(chB, :);
%     %             
%                 noRipFRa = sum(noRipMask(round(uTa))) / sum(noRipMask) * 1e3; %Hz
%                 noRipFRb = sum(noRipMask(round(uTb))) / sum(noRipMask)* 1e3; %Hz
%     % 
%                 uTnoRa = uTa(noRipMask(round(uTa)));
%                 uTnoRb = uTb(noRipMask(round(uTb)));
%     % 
%                 coF = 0;
%                 for ii = 1:length(uTnoRa)
%                     uT = uTnoRa(ii);
%                     if any(abs(uT - uTnoRb) <= coFireThresh); coF = coF + 1; end
%                 end
%     
%                 coFire = coF / sum(noRipMask) * 1e3;
%                 coFireChanceAna = 2 * (noRipFRa * noRipFRb) * coFireThresh/1e3;
%                 coFireChance = computeChanceCoFire(uTnoRa, uTnoRb, noRipMask, coFireThresh, 10);

%                 interactionType{uA,chBi} = sprintf('%s-%s', typeA, typeB);
            end
        end
    end
    
end
exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/CoFire';
save(fullfile(exportDir, [subject, '_', state, '_coFire-v2.mat']), 'coFireProbCoR', ...
                'coFireProbCoR_Shuff' , ...
                'coFireProbCoR_Ana', ...
                'coRippleFRa', ...
                'coRippleFRb', ...
                'coRippleDur', '-v7.3')

%%
% close all


exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/CoFire';

types = interactionType(:);
keep = cellfun(@(X) ~isempty(X), types);
types = types(keep);
types = unique(types);

figure;
for ii = 1:length(types)
    % boxplot(coFireProbCoR_2(coFireProbCoR_2 > 0), 'BoxStyle','filled');
    subplot(2,2,ii)
    ind = strcmp(interactionType,types(ii));
    ind = ind & coFireProbCoR > 0 & coFireProbCoR_Shuff > 0;
    bar([mean(coFireProbCoR(ind)./coFireProbCoR_Shuff(ind)) mean(coFireProbCoR(ind)./coFireProbCoR_Ana(ind))]); hold on;
    title(types(ii))
    hline(1)
end
%    h=findobj(gca,'tag','Outliers'); % Get handles for outlier lines.

filename = sprintf('%s_%s_coFire_10ms.mat',subject, state);
save(fullfile(exportDir, filename), '-v7.3')

%%

% coFireThreshAll = [-50:50];
unitChannels = unique(cell2mat(units(:,1)));
% CCG_coRip = cell(length(units),length(units));
% CCG_noRip = cell(length(units),length(units));
coRip_phaseA_all = cell(length(unitChannels),length(unitChannels));
coRip_phaseB_all = cell(length(unitChannels),length(unitChannels));
coRipPLV = nan(length(unitChannels),length(unitChannels));
noRipPLV = nan(length(unitChannels),length(unitChannels));

interactionType = cell(length(units),length(units));
% coRipPLV = [];
for chAi = 1:length(unitChannels)
    for chBi = 1:length(unitChannels)
        
        
        chA = unitChannels(chAi); chB = unitChannels(chBi);

        clc
        fprintf('computing PLV \n %i %i\n\n\n\n', chA, chB)
        
        if chA ~= chB 
%             uTa = units{uA,2} * 1e3;
%             uTb = units{uB,2} * 1e3;
%             
            coRipMask = rippMask(chA, :) & rippMask(chB, :);
            noRipMask = ~rippMask(chA,:) & ~rippMask(chB,:) & StageMask_1kHz;

%             noRipMask = ~rippMask(chA, :) & ~rippMask(chB, :) & StageMask_1kHz;
%             
% 
%             uTcoRa = uTa(coRipMask(round(uTa)));
% 
%             coRipMaskbnds = mask2bounds(coRipMask);
%             coRipMaskExtend = bounds2mask(coRipMaskbnds, length(rippMask), max(coFireThreshAll)); %pad by 25ms
%             uTcoRb = uTb(coRipMaskExtend(round(uTb)));
% 
%             uTnoRa = uTa(noRipMask(round(uTa)));
%             uTnoRb = uTb(noRipMask(round(uTb)));

%             sNo = uTnoRa - uTnoRb';
            
%             sNo = nan(1,length(uTnoRa)); %cross correlogram array
%             win = 100;
%             for s = 1:length(uTnoRa)
%                 spk = uTnoRa(s);
%                 if spk > win 
%                     ii = uTnoRb >= (spk-win) & uTnoRb < (spk+win);
%                     times = uTnoRb(ii) - spk;
%                     if size(times,1) > 1; times = times'; end
%                     if ~isempty(times)
% 
%                         sNo = [sNo times];
%                     end
%                 end
%             end
% 
%             sNo(isnan(sNo)) = [];
            
            
            
%             sCo = uTcoRa - uTcoRb';
% 
%             CCG_coRip{uA, uB} = sCo;
%             CCG_noRip{uA, uB} = sNo;

%             ii = arrayfun(@(X,Y) any((X > uTcoRa & Y < uTcoRa)) | any((X > uTcoRb & Y < uTcoRb)), coRipMask(:,2), coRipMask(:,1));
%             coRipMaskSpike = coRipMask(ii,:);
%             coRipMaskSpike = bounds2mask(coRipMaskSpike, length(rippMask));
%             phaseA = RBphaseAll(chA,coRipMask);
%             phaseB = RBphaseAll(chB,coRipMask);
% 
%             plvAB = PLV(phaseA,phaseB);
            
            
            
            phaseA = RBphaseAll(chA,coRipMask);
            phaseB = RBphaseAll(chB,coRipMask);
            coRip_phaseA_all{chA,chB} = phaseA;
            coRip_phaseB_all{chA,chB} = phaseB;
            plvAB = PLV(phaseA,phaseB);
            coRipPLV(chA,chB) = plvAB;

            
            phaseA = RBphaseAll(chA,noRipMask);
            phaseB = RBphaseAll(chB,noRipMask);

            plvAB = PLV(phaseA,phaseB);
            noRipPLV(chA,chB) = plvAB;
            
            
%             interactionType{chA,chB} = sprintf('%s-%s', typeA, typeB);

%             if ~isempty(uTcoRa) && ~isempty(uTcoRb)
%                 coRipFRa = sum(coRipMaskExtend(round(uTa))) / sum(coRipMaskExtend) * 1e3; %Hz
%                 coRipFRb = sum(coRipMaskExtend(round(uTb))) / sum(coRipMaskExtend) * 1e3; %Hz
%     
%                parfor thr = 1:length(coFireThreshAll)
%                     coFireThresh = coFireThreshAll(thr);
%                     
%                     s = uTcoRa - uTcoRb';
% %                     s = abs(s);
% 
%                     coFire = sum(s(:) <= coFireThresh & s(:) > (coFireThresh - 1));
%                     
%                     
%                     
% 
%                     CCG_All(uA,uB,thr) = coFire ;%/ coFireChanceAna;
%                     coRip_All(uA,uB,thr) = sum(coRipMaskExtend) * 1e3 ;%/ coFireChanceAna;
%         
%     %                 noRipMask = ~rippMask(chA, :) & ~rippMask(chB, :);
%     %     %             
%     %                 noRipFRa = sum(noRipMask(round(uTa))) / sum(noRipMask) * 1e3; %Hz
%     %                 noRipFRb = sum(noRipMask(round(uTb))) / sum(noRipMask)* 1e3; %Hz
%     %     % 
%     %                 uTnoRa = uTa(noRipMask(round(uTa)));
%     %                 uTnoRb = uTb(noRipMask(round(uTb)));
%     %     % 
%     %                 coF = 0;
%     %                 for ii = 1:length(uTnoRa)
%     %                     uT = uTnoRa(ii);
%     %                     if any(abs(uT - uTnoRb) <= coFireThresh); coF = coF + 1; end
%     %                 end
%     %     
%     %                 coFire = coF / sum(noRipMask) * 1e3;
%     %                 coFireChanceAna = 2 * (noRipFRa * noRipFRb) * coFireThresh/1e3;
%     %                 coFireChance = computeChanceCoFire(uTnoRa, uTnoRb, noRipMask, coFireThresh, 10);
%     
%                 end
%                 interactionType{uA,uB} = sprintf('%s-%s', typeA, typeB);
%                 coRipPLV(uA,uB) = plvAB;
%             end

            
        end
    end
%     uA
end
%%
% close all
% figure;
% for thr = 1:10
%     coFRatio = coFireProbCoR_Ana_All(:,:,thr);
%     coFRatio = coFRatio(coFRatio > 0);
%     plot(coFireThreshAll(thr),   mean(coFRatio), 'o'); hold on;
% end
 
load('T11_CCGs_NREM.mat', 'types', 'interactionType', 'coRipPLV', 'coRip_All', 'CCG_All', 'PLVbins', 'coFireProbCoR', 'coFireThreshAll')
% CCG_All(105:end,:,:) = [];
% coRip_All(105:end,:,:) = [];

for ii = 1:length(types)
    indType = strcmp(interactionType,types(ii));
%     indType = indType & coFireProbCoR > 0;
    
    yy = [];
    yySEM = [];
    %     for thr = 1:(length(coFireThreshAll))
    %         coFRatio = coFireProbCoR_Ana_All(:,:,thr);
    %         coFRatio = coFRatio(ind);
    %         coFRatio = coFRatio(coFRatio > 0);
    %         yy = [yy,   sum(coFRatio, 'omitnan')];
%         yySEM = [yySEM, std(coFRatio)/sqrt(length(coFRatio))];
%     end

    
    PLVbins = quantile(coRipPLV(indType), 100); %:0.01:0.9;
    diffPLV = [PLVbins(1) diff(PLVbins)];
    PLV_CCG = zeros(length(PLVbins), size(CCG_All,3));
    for b = 1:length(PLVbins)
        ind = coRipPLV > (PLVbins(b)-diffPLV(b)) & coRipPLV < PLVbins(b);
        CCGtemp = reshape(CCG_All, [size(CCG_All,1)*size(CCG_All,2) length(coFireThreshAll)]);
        timeTemp = reshape(coRip_All, [size(CCG_All,1)*size(CCG_All,2) length(coFireThreshAll)]);

        ind = reshape(ind, [size(coRipPLV,1)*size(coRipPLV,2)  1]);
        indType = reshape(indType, [size(coRipPLV,1)*size(coRipPLV,2)  1]);
        if sum(ind & indType) > 0
            PLV_CCG(b,:) = sum(CCGtemp(ind & indType,:)) / sum(timeTemp(ind & indType,1));
        end

%         sum(ind & indType)
    end

   

    tmp = coFireProbCoR(indType);
    tmpPLV = coRipPLV(indType);

    figure;

    % boxplot(coFireProbCoR_2(coFireProbCoR_2 > 0), 'BoxStyle','filled');
    subplot(2,2,ii)
%     yy = detrend(smoothdata(yy/max(yy), 'gaussian', 5));
    yySEM = yySEM/max(yy);
    yy = smoothdata(yy/max(yy), 'gaussian', 10);
%     yy = yy/max(yy);
%     xx = pwelch(yy/max(yy),1000);
    
%     plot(coFireThreshAll,yy); hold on;
%     boundedline(coFireThreshAll,yy,yySEM); hold on;
%     imagesc(PLV_CCG); hold on; colorbar;
%     imagesc(imgaussfilt(PLV_CCG(:,:),1)); hold on; colorbar;
    keep = sum(CCGtemp, 2, 'omitnan') > 0 & indType;
    mu = smoothdata(mean(CCGtemp(keep,:), 'omitnan'),'gaussian', 5);
    sigma = smoothdata(std(CCGtemp(keep,:), 'omitnan'),'gaussian', 5)/sqrt(sum(~isnan(CCGtemp(:,1))));
    boundedline(-50:50, mu, sigma);
    vline([-33:11:33])
%     dscatter(tmpPLV, tmp); hold on;
%     pl = plot(PLV_CCG');
%     cmap = colormap(autumn(length(PLVbins)))
%     for iP = 1:length(pl); pl(iP).Color = cmap(iP,:); end
%     colormap('autumn')
    box off
% axis off
% 
%     figure(2);
%     subplot(2,2,ii)
%     plot(PLVbins)



    title(types(ii))
end
%%
% coRipMask = zeros(96,96,length(rippMask));
% noRipMask = zeros(96,96,length(rippMask));
%                 for chA = 1:96
%     for chB = 1:96
%         coRipMask(chA,chB,:) = rippMask(chA, :) & rippMask(chB, :);
%         noRipMask(chA,chB,:) = ~rippMask(chA, :) & ~rippMask(chB, :) & StageMask_1kHz';
% 
%     end
% 
%     chA
% end
%% co fire during co ripple periods

coRipSpike = cell(length(units),length(units));
coRipDur = cell(length(units),length(units));
noRipSpike = cell(length(units),length(units));
interactionType = cell(length(units),length(units));
coFireProbNoR = cell(length(units),length(units));

for uAi = 1:length(units)
    for uBi = 1:length(units)
%         uA = unitsA(uAi); uB = unitsB(uBi);
        uA = uAi; uB = uBi;
        chA = units{uA,1}; chB = units{uB,1};
        typeA = units{uA,3}; typeB = units{uB,3};
        
        if chA ~= chB && any(strcmp(typeA,{'pyr', 'int'})) && any(strcmp(typeB,{'pyr', 'int'}))
            uTa = units{uA,2} * 1e3;
            uTb = units{uB,2} * 1e3;

            coRipMask = rippMask(chA, :) & rippMask(chB, :);
            noRipMask = ~rippMask(chA, :) & ~rippMask(chB, :) & StageMask_1kHz;

            

            uTcoRa = uTa(coRipMask(round(uTa)));
            uTcoRb = uTb(coRipMask(round(uTb)));

            uTnoRa = uTa(noRipMask(round(uTa)));
            uTnoRb = uTb(noRipMask(round(uTb)));

            controlMask = zeros(size(coRipMask));

            if ~isempty(uTcoRa) && ~isempty(uTcoRb)
                
                
                bnds = mask2bounds(coRipMask);
                dur = bnds(:,2) - bnds(:,1);
                bnds(dur < 25, :) = [];
                dur(dur < 25) = [];
                
                coSpikeCount = nan(2,length(bnds));
                noSpikeCount = nan(2,length(bnds));
                bndsBaselineAll = nan(size(bnds)); 
                for ibnd  = 1:length(bnds)
                    aSpikesCo = sum(uTcoRa >= bnds(ibnd,1)   &  uTcoRa <= bnds(ibnd,2));
                    bSpikesCo = sum(uTcoRb >= bnds(ibnd,1)   &  uTcoRb <= bnds(ibnd,2));

                    shift = 2e3;
                    bndsBaseline = bnds(ibnd,:) - shift;
                    if any(bndsBaseline < 0); continue; end

                    loopCount = 0;
                    while ~all(noRipMask(bndsBaseline(1):bndsBaseline(2))) && loopCount < 1000  
                        bndsBaseline = bndsBaseline - 100;
                        loopCount = loopCount + 1;
                        if any(bndsBaseline < 1)
                            loopCount = 1000;
                            break
                        end
                    end

                    if loopCount > 1000; continue; end

                    aSpikesNo = sum(uTnoRa >= bndsBaseline(1)   &  uTnoRa <= bndsBaseline(2));
                    bSpikesNo = sum(uTnoRb >= bndsBaseline(1)   &  uTnoRb <= bndsBaseline(2));

                    coSpikeCount(:, ibnd) = [aSpikesCo, bSpikesCo];
                    
                    noSpikeCount(:, ibnd) = [aSpikesNo bSpikesNo];
                    bndsBaselineAll(ibnd, :) = bndsBaseline;

                end

               controlMask = bounds2mask(bndsBaselineAll, length(coRipMask));
               uTcontrola = uTa(controlMask(round(uTa)));
               uTcontrolb = uTb(controlMask(round(uTb)));
    
                s = uTcontrola - uTcontrolb';
                N = histcounts(s, -50.5:50.5);

           
                coRipSpike{uA,uB} = coSpikeCount;
                coRipDur{uA,uB} = dur;
                noRipSpike{uA,uB} = noSpikeCount;
                coFireProbNoR{uA,uB} = N;

                interactionType{uA,uB} = sprintf('%s-%s', typeA, typeB);
            end
        end
    end
    uA
end

nCoSpikeCoAll = [];
nCoSpikeNoAll = [];

distanceAll = [];
c = 1;

nCoSpikeNo = [];
nCoSpikeCo = [];
%     dur = 0;
for uAi = 1:length(units)
    for uBi = 1:length(units)
%         uA = units(uAi); uB = units(uBi);
        uA = uAi; uB = uBi;
        chA = units{uA,1};
        chB = units{uB,1};
        d = findContactDistance(chMap, chA,chB);
        
            coRip = coRipSpike{uA,uB};
            noRip = noRipSpike{uA,uB};
            dur = sum(coRipDur{uA,uB}) / 1e3;

            if ~isempty(coRip) 
                sp = coRip(:,coRip(1,:) > 0 & coRip(2,:) > 0);
                nCoSpikeCo = size(sp,2); %/dur; % / (dur);
                
                sp = noRip(:,noRip(1,:) > 0 & noRip(2,:) > 0);
                nCoSpikeNo = size(sp,2); %/dur; % / (dur);
%                     dur = dur + (sum(coRipDur{uA,uB}) / 1e3);

            nCoSpikeCoAll = [nCoSpikeCoAll, nCoSpikeCo];
            nCoSpikeNoAll = [nCoSpikeNoAll, nCoSpikeNo];
            distanceAll = [distanceAll, d];
            end
        
            

    end
end

%% Calculate Distance 
home_dir = '/space/seh8/2/mdeh1-3/halgdev/projects/cgonzalez/Units';
chmap_file = sprintf('%s/mg49_and_mg29_neuroport_channel_map.csv',home_dir);
chMap = dlmread(chmap_file);

distances = nan(length(units), length(units));
for uA = 1:length(units)
    for chBi = 1:length(units)
        chAi = units{uA,1};
        chB = units{chBi,1};
        distances(uA,chBi) = findContactDistance(chMap, chAi,chB);

    end
end















%%


nCoSpikeCoDist = [];
nCoSpikeNoDist = [];

nCoSpikeCoDistSEM = [];
nCoSpikeNoDistSEM = [];


bins = [0:500:4500]; %+(8*400);
bins = bins/400;
dW = diff(bins);
bins(end) = [];

c = 1;

for W = bins
    nCoSpikeNo = [];
    nCoSpikeCo = [];
%     dur = 0;
    for uAi = 1:length(unitsA)
        for uBi = 1:length(unitsA)
            uA = unitsA(uAi); chBi = unitsA(uBi);
            chAi = units{uA,1};
            chB = units{chBi,1};
            d = findContactDistance(chMap, chAi,chB);
            
            if d > W && d < W + dW(c)
                coRip = coRipSpike{uA,chBi};
                noRip = noRipSpike{uA,chBi};
                dur = sum(coRipDur{uA,chBi}) / 1e3;
    
                if ~isempty(coRip) 
                    sp = coRip(:,coRip(1,:) > 0 & coRip(2,:) > 0);
                    nCoSpikeCo = [nCoSpikeCo, size(sp,2)/dur]; % / (dur);
                    
                    sp = noRip(:,noRip(1,:) > 0 & noRip(2,:) > 0);
                    nCoSpikeNo = [nCoSpikeNo, size(sp,2)/dur]; % / (dur);
%                     dur = dur + (sum(coRipDur{uA,uB}) / 1e3);
                end
            end
                
    
        end
    end

    nCoSpikeCoDist = [nCoSpikeCoDist, mean(nCoSpikeCo)];
    nCoSpikeNoDist = [nCoSpikeNoDist, mean(nCoSpikeNo)];

    nCoSpikeCoDistSEM = [nCoSpikeCoDistSEM, computeSEM(nCoSpikeCo, 1)];
    nCoSpikeNoDistSEM = [nCoSpikeNoDistSEM, computeSEM(nCoSpikeNo,1)];

    c = c + 1
end
%%
figure('Position', [995 908 306 413]);
subjects = {'T11', 'T11-dual', 'MG29', 'MG63'};
state = 'NREM';
clr = [237/255 110/255 87/255; %T11
       255/255 66/255 161/255; %T11-dual
       108/255 227/255 208/255; %MG29
       129/255 213/255 83/255]; %MG63

clr2 = brewermap(10,'Paired');

p = polyfit(1:length(nCoSpikeCoDist), nCoSpikeCoDist,2); %quadratic fit
x1 = linspace(1,length(nCoSpikeCoDist));
y1 = polyval(p,x1);
l1 = plot(x1,y1); hold on;
l1.Color =  clr2(8,:);
b1 = errorbar(1:length(nCoSpikeCoDist), nCoSpikeCoDist, nCoSpikeCoDistSEM, 'o'); hold on;
b1.MarkerFaceColor = clr(4,:);
b1.Color = clr2(8,:);
b1.LineWidth = 1.5;




p = polyfit(1:length(nCoSpikeNoDist), nCoSpikeNoDist,2); %quadratic fit
x1 = linspace(1,length(nCoSpikeNoDist));
y1 = polyval(p,x1);
l2 = plot(x1,y1); hold on;
l2.Color =  clr2(1,:);
b2 = errorbar(1:length(nCoSpikeNoDist), nCoSpikeNoDist, nCoSpikeNoDistSEM, 'o'); hold on;
b2.MarkerFaceColor = clr(4,:);
b2.Color = clr2(1,:);
b2.LineWidth = 1.5;


ylabel('CoFire Event Rate [Hz]')
% legend('CoRipples','Control')

% yyaxis right
% pl = plot(nCoSpikeCoDist./nCoSpikeNoDist);
% pl.Color = [180/255 160/255 101/255];
% pl.LineWidth = 2;
% title(state)
fig = gcf;
fig.Color = 'w';
box off

xlim([0 10])
ylim([0 1.2])

%% CoRipple Co Spike avalance
if exist('rippMedial','var'); rippleStats = rippMedial.rippleStats; end
RippleSpikeProb = zeros(1,96);
RippleSpikeMask = zeros(size(rippMask));
unitChannels = units(cellfun(@(X) strcmp(X,'pyr') | strcmp(X,'int'), units(:,3)),1); 
unitChannels = unique(cell2mat(unitChannels))';
for iCh = unitChannels
    if iCh <= 96
        iStart = rippleStats.window{iCh}(:,1);
        iEnd = rippleStats.window{iCh}(:,2);
    else
        iStart = rippLateral.rippleStats.window{iCh-96}(:,1);
        iEnd = rippLateral.rippleStats.window{iCh-96}(:,2);
    end
    chSpikes = spikeMask(iCh,:);
       
    c = 0;
    for r = 1:length(iStart)
        if any(chSpikes(iStart(r):iEnd(r)))
           c = c + 1; 
           RippleSpikeMask(iCh,iStart(r): iEnd(r)) = 1;
        end
    end
        
%     RippleSpikeProb(iCh) = c / sum(rippMask(iCh,:)); % control by duration     
    RippleSpikeProb(iCh) = c / length(iStart); % control by nuber of ripples   
            
    
end

nCoRipSpike = sum(RippleSpikeMask);

CoupledAll = nan(1,1e6);
unCoupledAll = nan(1,1e6);
nAll = nan(1,1e6);
ii = 1;
for n = 2:20
    BoundsCoRip = mask2bounds(nCoRipSpike == n);
    
    if size(BoundsCoRip,1) > 5e2
        xx = randsample(size(BoundsCoRip,1), 5e2)';
    else
        xx = 1:size(BoundsCoRip,1);
    end

    
    for iB = xx
       
       ch = find(RippleSpikeMask(:,BoundsCoRip(iB,1)));
       unCoupled = prod(RippleSpikeProb(ch));
       
       boundCoRipple = mask2bounds(sum(rippMask(ch,:)) == n);
       boundCoRippleSpike = mask2bounds(sum(RippleSpikeMask(ch,:)) == n);
       
       unCoupledAll(ii) = unCoupled;
       CoupledAll(ii) = size(boundCoRippleSpike,1) / size(boundCoRipple,1);
       nAll(ii) =  n;

       ii = ii + 1;

    end

    n
    
    
end




%%
close all
RatioBinned = [];
RatioSEM = [];

for n = 2:20
    ii = find(nAll == n);
    ratio = log10(CoupledAll(ii)./unCoupledAll(ii));
    
    RatioBinned = [RatioBinned, mean(ratio)];
    RatioSEM = [RatioSEM, std(ratio)/sqrt(length(ii))];

    
    
end

h4  =  figure('Position',[629 426 507 279]);
boundedline(2:20,RatioBinned,RatioSEM); hold on;
% save(fullfile(fdir,'MG29_CoRipleCoSpike.mat'), 'coSpikeMatrixCoupled', 'coSpikeMatrixUncoupled', '-v7.3')
xlabel('# co Occuring ripples with spikes')
ylabel('log_{10}(observed prob / unCoupled prob)')

ax = gca;
ytick = ax.YTick;
ytick = num2cell(ytick);
ax.YTickLabel = ytick;

fig = gcf;
fig.Color = 'w';





