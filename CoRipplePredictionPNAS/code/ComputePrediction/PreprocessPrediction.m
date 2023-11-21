clear
close all
clc
addpath(genpath('../../code'))

% parpool(12)
%%
subject = 'E1';
state = 'wake';

outputDir = sprintf('./out/%s', subject);
if ~isfolder(outputDir); mkdir(outputDir); end

units = load(); %see readMe for instructions on how to format unit matrix

load() %load rippleband phase data

matExportFolder = ''; %folder for output of ripple detection.
load(fullfile(matExportFolder,filename));

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
computeSingleRipPredict = true; 
% true - run anlysis when only A unit has a ripple as control 
% false - run anlysis when no unit has a ripple as control 

CCGfolder = sprintf('./out/%s/CCG/%s', subject, state); %folder where output from computeCoRippleCCG.m


win = 200;
nIter = 1;
smoothingWindow = 25;   
minSpikes = 15;

nullPredictAll = [];
rippPredictAll = [];

nullPhaseAll = [];
rippPhaseAll = [];

rippFiringRate = [];
nullFiringRate = [];

rippBins = [];
nullBins = [];

distances = [];
coRippleRate = [];

MannUAll = [];
diffAll = [];   

interactionType = {};

rippStats = [];
nullStats = [];

for a = 1:length(unitsAll)
    clear ripp; clear null;
    ua = unitsAll(a);
    
    chA = find(rippleChannels == units{ua,1});
    typeA = units{ua,3};
    map = chMap == chA;
    [Xaii, Yaii] = find(map);

    rMaska = rippMask(chA,:);

    fname = sprintf('%s_CCGs_%s_unit_%i.mat', subject, state, a);
    
    
    clc
    fprintf('... loading %s ...\n', fname)
    load(fullfile(CCGfolder,fname))
    fprintf('... done ...\n')

    if computeSingleRipPredict
        null = rippA;
    end

    for b = 1:length(unitsAll)

        if a == b
            continue
        end

        ub = unitsAll(b);
    
        chB = find(rippleChannels == units{ub,1});
        typeB = units{ub,3};

        rMaskb = rippMask(chB,:);

        nControl = sum(~isnan(null.History{a,b}(:,1)));
        nRipple  = sum(~isnan(ripp.History{a,b}(:,1)));

        if chA == chB || nRipple == 0
            continue
        end
        
        

        
        rippHistory = ripp.History{a,b}(:,end-(win-1):end);
        nullHistoryAll = null.History{a,b}(:,end-(win-1):end);

        nRippleSpikes = sum(rippHistory(:), 'omitnan');
        fprintf('... chan A %.02f %%  ... chan B %.02f %%  ...\n', a/length(unitsAll)*100, b/length(unitsAll)*100)
        fprintf('... nRippl %.02i     ... nSpike %.02i  ...\n', nRipple, nRippleSpikes)
        
        if (nRippleSpikes > minSpikes) && nRipple > 2

            for iter = 1:nIter

                nLoop = 0;

                if (sum(nullHistoryAll(:), 'omitnan') > sum(rippHistory(:), 'omitnan')) && nRipple < nControl
                    Y = randsample(nControl, nRipple); %randomly sample from control distribution;
                    nullHistory = null.History{a,b}(Y,end-(win-1):end);
                    while sum(nullHistory(:), 'omitnan') < sum(rippHistory(:), 'omitnan') && nLoop < 1000    
                        newY = randsample(nControl, 1);
                        while ismember(newY,Y)
                            newY = randsample(nControl, 1);
                        end
    
                        Y = [Y; newY];
                        rippHistory = ripp.History{a,b}(:,end-(win-1):end);
                        nullHistory = null.History{a,b}(Y,end-(win-1):end);
    
                        nLoop = nLoop + 1;
                                        
                    end
                else
                    Y = 1:nControl;
                end

                fprintf('nLoop = %i\n', nLoop)

                rippHistory = ripp.History{a,b}(:,end-(win-1):end);
                nullHistory = null.History{a,b}(Y,end-(win-1):end);

                rippNan = isnan(rippHistory(:,1));
                nullNan = isnan(nullHistory(:,1));
                rippHistory(rippNan,:) = [];
                nullHistory(nullNan,:) = [];

                
                
                nullHistoryFilter = null.History{a,b}(:,2:end);
                rippHistoryFilter = ripp.History{a,b}(:,2:end);

                rippNan = isnan(rippHistoryFilter(:,1));
                nullNan = isnan(nullHistoryFilter(:,1));
                rippHistoryFilter(rippNan,:) = [];
                nullHistoryFilter(nullNan,:) = [];


                APar = ripp.History{a,b}(~rippNan,1);
                APan = null.History{a,b}(:,1);
                APan(nullNan) = [];

                nullPredict = [];
                nullPhaseA = [];
                nullPhaseB = [];
                nullTimes = [];

                [nullTrials, nullTimes] = find(nullHistoryFilter);

                [rippTrials, rippTimes] = find(rippHistoryFilter);

                rippStats.times{a,b} = rippTimes;
                rippStats.trials{a,b} = rippTrials;
                rippStats.APar{a,b} = APar;
                rippStats.nTrials(a,b) = size(rippHistoryFilter,1);
                rippStats.unitsA = unitsA;
                rippStats.unitsB = unitsB;

                nullStats.times{a,b} = nullTimes;
                nullStats.trials{a,b} = nullTrials;
                nullStats.APan{a,b} = APan;
                nullStats.nTrials(a,b) = size(nullHistoryFilter,1);
                nullStats.unitsA = unitsA;
                nullStats.unitsB = unitsB;

                nullHistory(isnan(nullHistory)) = 0;
                nullHistoryFilter(isnan(nullHistoryFilter)) = 0;

                for iH = 1:size(nullHistory, 1)

                   nullHistorySmooth = smoothdata(sum(nullHistory(1:end ~= iH,:)), 'gaussian',smoothingWindow);
                   nullHistoryFilterSmooth = smoothdata(sum(nullHistoryFilter(1:end ~= iH,:)), 'gaussian',smoothingWindow);
                   
                   nullFilter =  constructFilter(nullHistorySmooth, nullHistoryFilterSmooth);
                   if length(nullFilter) == 1
                       nullFilter = zeros(1,win);
                   else
                       nullFilter = nullFilter(end-(win-1):end);
                   end
                   
                   iPhaseA = RBphaseAll(chA,round(APan(iH))-win+1:round(APan(iH)));
                   iPhaseB = RBphaseAll(chB,round(APan(iH))-win+1:round(APan(iH)));
                   if ~any(isnan(nullHistory(iH,:)))
                       filterOut = nullFilter(logical(nullHistory(iH,end-(win-1):end)));
                       phaseOutB = iPhaseB(logical(nullHistory(iH,end-(win-1):end)));
                       phaseOutA = iPhaseA(logical(nullHistory(iH,end-(win-1):end)));
                       times = find(nullHistory(iH,end-(win-1):end));

                       
                   
                       if isempty(filterOut); filterOut = nan; phaseOutB = nan; times = nan; end
                       if any(~isfinite(filterOut)); filterOut = nan; phaseOutB = nan; times = nan; end
                   else
                       filterOut = nan; phaseOutB = nan;
                   end
                   
                   nullPredict = [nullPredict, filterOut];
                   nullPhaseB = [nullPhaseB, phaseOutB];
                   nullPhaseA = [nullPhaseA, phaseOutA];
                   nullTimes = [nullTimes, times];

                end
               
                rippPredict = [];
                rippPhaseA = [];
                rippPhaseB = [];
                rippTimes = [];
                phaseDelays = [];
                
                rippHistory(isnan(rippHistory)) = 0;
                rippHistoryFilter(isnan(rippHistoryFilter)) = 0;
                for iH = 1:size(rippHistory, 1)

                   rippHistorySmooth = smoothdata(sum(rippHistory(1:end ~= iH,:)), 'gaussian',smoothingWindow);
                   rippHistoryFilterSmooth = smoothdata(sum(rippHistoryFilter(1:end ~= iH,:)), 'gaussian',smoothingWindow);

                   rippFilter =  constructFilter(rippHistorySmooth, rippHistoryFilterSmooth);
                   rippFilter = rippFilter(end-(win-1):end);
                     
                
                   iPhaseA = RBphaseAll(chA,round(APar(iH))-win+1:round(APar(iH)));
                   iPhaseB = RBphaseAll(chB,round(APar(iH))-win+1:round(APar(iH)));

                   if ~any(isnan(rippHistory(iH,:)))
                       phaseOutB = iPhaseB(logical(rippHistory(iH,end-(win-1):end)));
                       phaseOutA = iPhaseA(logical(rippHistory(iH,end-(win-1):end)));
                       filterOut = rippFilter(logical(rippHistory(iH,end-(win-1):end)));
                       times = find(rippHistory(iH,end-(win-1):end));

                       pDelay = phaseOutB - phaseOutA;

                       if isempty(filterOut); filterOut = nan; phaseOutB = nan; times = nan; pDelay = nan; end
                       if any(~isfinite(filterOut)); filterOut = nan; phaseOutB = nan; times = nan; end

                   else
                       filterOut = nan; phaseOutB = nan; 
                   end
                   
                   rippPredict = [rippPredict, filterOut];
                   rippPhaseB = [rippPhaseB, phaseOutB];
                   rippPhaseA = [rippPhaseA, phaseOutA];
                   rippTimes = [rippTimes, times];
                   phaseDelays = [phaseDelays, pDelay];

                end


                rippStats.FilterOut{a,b} = rippPredict;
                nullStats.FilterOut{a,b} = nullPredict;

                rippStats.nBaseline{a,b} = sum(rippHistoryFilter(:));
                nullStats.nBaseline{a,b} = sum(nullHistoryFilter(:));

                
                rippStats.PhaseA{a,b} = rippPhaseA;
                nullStats.PhaseA{a,b} = nullPhaseA;
                rippStats.PhaseB{a,b} = rippPhaseB;
                nullStats.PhaseB{a,b} = nullPhaseB;
                rippStats.PhaseDelay{a,b} = phaseDelays;
    
                phases = RBphaseAll(chA, round(APar));
                rippStats.PLV_A(a,b) = PLV(phases, zeros(1,length(phases)));
                phases = RBphaseAll(chA, round(APan));
                nullStats.PLV_A(a,b) = PLV(phases, zeros(1,length(phases)));

                rippStats.times{a,b} = rippTimes;
                nullStats.times{a,b} = nullTimes;

                rippStats.Filter{a,b} = rippFilter;
                nullStats.Filter{a,b} = nullFilter;
                
                [~, h] = signrank(rippPredict);
                rippStats.Sig(a,b) = h;
                [~, h] = signrank(nullPredict);
                nullStats.Sig(a,b) = h;

                nullPredictAll = [nullPredictAll, median(nullPredict, 'omitnan')];
                rippPredictAll = [rippPredictAll, median(rippPredict, 'omitnan')];
                
%                 
                nControl = length(Y);
                nRipple  = size(ripp.History{a,b},1);

                rippHistory = ripp.History{a,b}(:,end-(win-1):end);
                nullHistory = null.History{a,b}(Y,end-(win-1):end);

                nControlNonZero = sum(sum(nullHistory,2) > 0);
                nRippleNonZero = sum(sum(rippHistory,2) > 0);

                
                rHz = sum(rippHistory(:)) / (win * nRippleNonZero) * rippleStats.fs;
                nHz = sum(nullHistory(:)) / (win * nControlNonZero) * rippleStats.fs;
                rippFiringRate = [rippFiringRate, rHz];
                nullFiringRate = [nullFiringRate, nHz];
               
                
                [~, times] = find(rippHistory);
                [N,~] = histcounts(times, 0:1:win);
                rippBins = [rippBins, N];

                [~, times] = find(nullHistory);
                [N,~] = histcounts(times, 0:1:win);
                nullBins = [nullBins, N];

                distances = [distances, findContactDistance(chMap, chA,chB)];
                
                bnds = mask2bounds(rMaska + rMaskb > 1);
                coRippleRate = [coRippleRate, size(bnds,1) / length(rippleStats.locs{chA})];

                [p,H,stats] = ranksum(rippPredict,nullPredict);
                

                MannUAll = [MannUAll, H];

                interaction = sprintf('%s->%s',typeB,typeA);
                rippStats.CellType{a,b} = interaction;
                interactionType{end+1} = interaction;

                
            end
        
        end


    end
    
     
    

end

[H,P,CI] = ttest(nullPredict2,rippPredict2)
mean(nullPredict2, 'omitnan')
mean(rippPredict2, 'omitnan')

if computeSingleRipPredict
    save(fullfile(outputDir,sprintf('Prediction_%s_%s-preprocessed-rippA.mat', subject, state)), 'nullStats', 'rippStats', '-v7.3')
else
    save(fullfile(outputDir,sprintf('Prediction_%s_%s-preprocessed.mat', subject, state)), 'nullStats', 'rippStats', '-v7.3')
end






