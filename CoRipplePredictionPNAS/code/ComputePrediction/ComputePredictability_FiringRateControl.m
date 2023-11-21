
clear
close all
clc

addpath(genpath('../../code'))

%%
% Load ouput from PreprocessPrediction.m
load(fullfile(outputDir,sprintf('Prediction_%s_%s-preprocessed.mat', subject, state)), 'nullStats', 'rippStats')

load() %load rippleband phase data

close all
subject = 'E1';
state = 'wake';
units = load(); %see readMe for instructions on how to format unit matrix

%loop through filter parameters

for win = 150 %0:50:600
%     win = 150;
    smoothingWindow = 5;
    smoothingWindowBin = 1;
    smoothingWindowExtend = 0;
    minSpikes = 25;

    unitsA = rippStats.unitsA;
    unitsB = rippStats.unitsB;
    unitsAll = sort(unique([unitsA;unitsB]));


    filterMode = 'differentBaseline'; %'differentBaseline'
    % Clear outputs
    
    rippPredictAll = [];
    nullPredictAll = [];
    distances = [];
    MannUAll = [];
    interactionType = {};
    rippStats.spearmanPval = [];
    rippStats.spearmanRho = [];
    
    rippStats.badAP = [];
    
    rippStats.spearmanRho = [];
    rippStats.spearmanPval = [];
    
    rippStats.timesWin = [];
    
    rippStats.predict = [];
    nullStats.predict = [];
    
    rippStats.FilterOut = [];
    nullStats.FilterOut = [];
    
    rippStats.Filter = [];
    nullStats.Filter = [];
    
    nullStats.Y = [];
    
    rippStats.FiringRate = [];
    nullStats.FiringRate = [];
    
    
    rippStats.Sig = [];
    nullStats.Sig = [];
    
    rippStats.CellType = [];
    
    rippStats.PhaseA = [];
    nullStats.PhaseA = [];
    rippStats.PhaseB = [];
    nullStats.PhaseB = [];
    
    rippStats.PLV_A = [];
    nullStats.PLV_A = [];
    
    
    rippStats.PLV_B = [];
    nullStats.PLV_B = [];
    
    rippStats.rippPhaseUnitA = [];
    
    rippStats.chA = [];
    rippStats.chB = [];
    
    % Compute Predictability Metrics
    h = waitbar(0,'computing predictability metrics ...');
    for a = 1:size(rippStats.APar,1)
        ua = unitsAll(a);
        typeA = units{ua,3};
    
        chA = units{ua,1};
        for b = 1:size(rippStats.APar,2)
            ub = unitsAll(b);
            typeB = units{ub,3};
            chB = units{ub,1};
           
            APar = rippStats.APar{a,b};

            filterOutRipp = [];
            timesRipp = [];
            rippPhaseA = [];
            rippPhaseB = [];
            rippPhaseUnitA = [];

            trials = rippStats.trials{a,b};

            timesAllRipp = 751 - rippStats.times{a,b}(ismember(rippStats.trials{a,b},trials));
            rippFiringRate = sum(timesAllRipp <= win) / (win/1000) / length(unique(trials(timesAllRipp <= win)));

            nSpikeRipp = sum(timesAllRipp <= win);
            timesAllNull = 751 - nullStats.times{a,b}(~ismember(nullStats.trials{a,b},badAP));
            trialsAllNull = nullStats.trials{a,b}(~ismember(nullStats.trials{a,b},badAP));
            nSpikeNull = sum(timesAllNull <= win);  
            nullFiringRate = 0;
            nLoop = 0;
            
            controlTrials = unique(nullStats.trials{a,b}(~ismember(nullStats.trials{a,b},badAP)));
            if nSpikeNull > nSpikeRipp 
                Y = randsample(controlTrials, 1); %randomly sample from control distribution;
                ii = ismember(trialsAllNull,Y);
                timesNull = timesAllNull(ii);
                nSpikeNull = sum(timesNull <= win);
                nullFiringRate = 0;
                add = 1;
                while nSpikeNull < nSpikeRipp
                    Y = randsample(controlTrials, 1+add); %randomly sample from control distribution;
                    ii = ismember(trialsAllNull,Y);
                    timesNull = timesAllNull(ii);
                    nSpikeNull = sum(timesNull <= win);
                    add = add + 1;

                    nullFiringRate = sum(timesAllNull(ii) <= win) / (win/1000) /length(unique(trialsAllNull(timesAllNull <= win & ii)));

                end

            else
                Y = 1:length(controlTrials);
                 
            end
            
            trialsRippUnique = unique(trials);
            trialsRipp = trials;
            sub = 0;
            while rippFiringRate > nullFiringRate && sub < length(trialsRippUnique)
                    Yripp = randsample(trialsRippUnique, length(trialsRippUnique)-sub); %randomly sample from ripple distribution;
                    
                    trialsRipp = trials(ismember(trials, Yripp));
                    timesRippShuff = timesAllRipp(ismember(trials, Yripp));

                    rippFiringRate = sum(timesRippShuff <= win) / (win/1000) / length(unique(trialsRipp(timesRippShuff <= win)));

                    sub = sub + 1;

            end

            trialsUnique = unique(trialsRipp);
            s = 0;
            for iH = 1:length(trialsUnique)
                if any(trials == trialsUnique(iH))
                    timeTrial = rippStats.times{a,b}(trials == trialsUnique(iH));                
                    timeTrial = 751 - timeTrial + 1;
                    timeTrial(timeTrial > win) = [];
                    s = s  + length(timeTrial);
    
                    times = rippStats.times{a,b}(ismember(trials, trialsUnique) & trials ~= trialsUnique(iH));
                    times = 751 - times;
    
                    switch filterMode
                        case 'commonBaseline' %needs fixing
                            [rippHistory, EDGES,BIN] = histcounts(times(times <= (win)), -0.5:smoothingWindowBin:(win+1), 'Normalization','count');
                            rippHistory = repelem(rippHistory,smoothingWindowBin);
                            
                            timesRippNull = [times; timesNull];
                            baselineFilter = histcounts(timesRippNull, win:win+smoothingWindowExtend, 'Normalization','count');
                            
                            rippHistorySmooth = smoothdata(rippHistory, 'gaussian',smoothingWindow);
                            baselineFilterSmooth = smoothdata(baselineFilter, 'gaussian',smoothingWindow);
    
                            baselineFilterSmooth = baselineFilterSmooth / (nTrialsNull + length(unique(trials))-1);
                            rippHistorySmooth = rippHistorySmooth / (length(unique(trials))-1);
    
                            rippFilter =  constructFilter(rippHistorySmooth, baselineFilterSmooth);
                        case 'differentBaseline'
                            [rippHistory, EDGES,BIN] = histcounts(times(times <= (win+smoothingWindowExtend)), -0.5:smoothingWindowBin:(win+smoothingWindowExtend+1), 'Normalization','count');
                            rippHistory = repelem(rippHistory,smoothingWindowBin);
                            rippHistorySmooth = zscore(smoothdata(rippHistory, 'gaussian',smoothingWindow));
                            rippFilter =  rippHistorySmooth(1:win+1);
                    end
                    
    
    
                    filterOutRipp = [filterOutRipp , rippFilter(timeTrial)];
                    
                    iPhaseA = fliplr(RBphaseAll(chA,round(APar(iH))-win+1:round(APar(iH)))); %flip to align with relative spike times
                    iPhaseB = fliplr(RBphaseAll(chB,round(APar(iH))-win+1:round(APar(iH))));
                    timesRipp = [timesRipp timeTrial'];
                    rippPhaseB = [rippPhaseB iPhaseB(timeTrial)];
                    rippPhaseA = [rippPhaseA iPhaseA(timeTrial)];
                    rippPhaseUnitA = [rippPhaseUnitA, repmat(iPhaseA(1), [1 length(timeTrial)])];
                    
                end
            end
            
            
    
            nRipple = length(unique(rippStats.trials{a,b}(~ismember(rippStats.trials{a,b},badAP))));
    
            rippPredict = mean(filterOutRipp, 'omitnan');
    
            
            nSpikeRipp = sum(timesAllRipp <= win);
            filterOutNull = [];
    
            if ~isnan(rippPredict) && nSpikeRipp > minSpikes
                APan = nullStats.APan{a,b};
                
    
                for iter = 1:nIter                    
            
                    if nSpikeNull > 0
                        APan = APan(Y);
                        
                        ii = ismember(trialsAllNull,Y);
                        
                        filterOutNullIter = [];
        
                        nullPhaseA = [];
                        nullPhaseB = [];
                        trials = trialsAllNull(ii);
                        timesAllNull = timesAllNull(ii);
                        for iH = 1:length(Y)
                            if any(trials == Y(iH))
                                timeTrial = timesAllNull(trials == Y(iH)) + 1;
                                timeTrial(timeTrial > win) = [];
                                
                                times = timesAllNull(trials ~= Y(iH));
                                times = times + 1;
                                
                       
                                [nullHistory, EDGES,BIN] = histcounts(times(times <= (win+smoothingWindowExtend)), -0.5:smoothingWindowBin:(win+smoothingWindowExtend+1), 'Normalization','count');
                                nullHistory = repelem(nullHistory,smoothingWindowBin);
                                nullHistorySmooth = zscore(smoothdata(nullHistory, 'gaussian',smoothingWindow));
                                nullFilter =  nullHistorySmooth(1:win+1);
                                
            
                
                                filterOutNullIter = [filterOutNullIter , nullFilter(timeTrial)];
                                
                                iPhaseA = fliplr(RBphaseAll(chA,round(APan(iH))-win+1:round(APan(iH)))); %flip to align with relative spike times
                                iPhaseB = fliplr(RBphaseAll(chB,round(APan(iH))-win+1:round(APan(iH))));
                                nullPhaseB = [nullPhaseB iPhaseB(timeTrial)];
                                nullPhaseA = [nullPhaseA iPhaseA(timeTrial)];
    
                                
                                
                            end
                        end
                        nullFiringRate = sum(timesAllNull <= win) / (win/1000) /length(unique(trials(timesAllNull <= win)));
                       
                        filterOutNull = [filterOutNull, filterOutNullIter];
            
                        
                    end
                    
                end
            end
            
            nullPredict = mean(filterOutNull, 'omitnan');
            if ~isempty(filterOutRipp) && ~isempty(filterOutNull) && all(isfinite([rippPredict, nullPredict])) && length(filterOutRipp) > minSpikes
    
                
                [H,P,CI,STATS] = ttest2(filterOutRipp,filterOutNull);
        
                MannUAll = [MannUAll, H];
                rippPredictAll = [rippPredictAll, rippPredict];
                nullPredictAll = [nullPredictAll, nullPredict];
                
                rippStats.timesWin{a,b} = timesRipp;
                rippStats.predict{a,b} = rippPredict;
                nullStats.predict{a,b} = nullPredict;
                
                rippStats.FilterOut{a,b} = filterOutRipp;
                nullStats.FilterOut{a,b} = filterOutNull;
    
                rippStats.Filter{a,b} = rippFilter;
                nullStats.Filter{a,b} = nullFilter;
    
                nullStats.Y{a,b} = Y;
        
                rippStats.FiringRate(a,b) = rippFiringRate;
                nullStats.FiringRate(a,b) = nullFiringRate;
                close all
    
        
                if (chA > 96 && chB <= 96) || (chA <= 96 && chB > 96)
                    distanceAdd = 8; %mm;
                else
                    distanceAdd = 0; %mm
                end

                rippStats.Sig(a,b) = ttest(filterOutRipp);
                nullStats.Sig(a,b) = ttest(filterOutNull);

                interaction = sprintf('%s->%s',typeB,typeA);
                rippStats.CellType{a,b} = interaction;
                interactionType{end+1} = interaction;
    
                rippStats.PhaseA{a,b} = rippPhaseA;
                nullStats.PhaseA{a,b} = nullPhaseA;
                rippStats.PhaseB{a,b} = rippPhaseB;
                nullStats.PhaseB{a,b} = nullPhaseB;
    
                phases = RBphaseAll(chA, round(APar));
                rippStats.PLV_A(a,b) = PLV(phases, zeros(1,length(phases)));
                phases = RBphaseAll(chA, round(APan));
                nullStats.PLV_A(a,b) = PLV(phases, zeros(1,length(phases)));
    
    
                phases = RBphaseAll(chB, round(APar));
                rippStats.PLV_B(a,b) = PLV(phases, zeros(1,length(phases)));
                phases = RBphaseAll(chB, round(APan));
                nullStats.PLV_B(a,b) = PLV(phases, zeros(1,length(phases)));

                rippStats.rippPhaseUnitA{a,b} = rippPhaseUnitA;
    
                rippStats.chA(a,b) = chA;
                rippStats.chB(a,b) = chB;
            end
        end
    
        waitbar(a/length(unitsAll),h)
    
    end
    
    close(h)
    
    [H,P,CI] = ttest(nullPredictAll,rippPredictAll);

    fprintf('filter window: %i ', win)
    fprintf('baseline window: %i \n', smoothingWindowExtend)
    fprintf('ripple vs null significant? %i \n', H);
    fprintf('ripple prediction %0.4f\n',mean(rippPredictAll, 'omitnan'));
    fprintf('null prediction %0.4f\n',mean(nullPredictAll, 'omitnan'));
    
    filename = sprintf('Prediction_%s_%s-processed-FRcontrol-%i_PreWin-%i_win.mat', subject, state, smoothingWindowExtend, win);
    save(filename, 'rippStats','nullStats', 'win','smoothingWindow','smoothingWindowBin','smoothingWindowExtend', 'minSpikes','state','-v7.3')
end





