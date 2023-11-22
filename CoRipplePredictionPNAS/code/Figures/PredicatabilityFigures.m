close all
clear
clc

addpath(genpath('../../code'))
%%
clr2 = brewermap(10,'Paired');  
subjects = {'B1', 'B1-dual', 'E1', 'E2'};
state = 'NREM';
clr = [237/255 110/255 87/255; %B1
       255/255 66/255 161/255; %B1-dual
       108/255 227/255 208/255; %E1
       129/255 213/255 83/255]; %E2

exportDirec = '../../figures';
if ~isfolder(exportDirec); mkdir(exportDirec); end
%% Predictability vs Distance (Fig. 4c)
% close all
for subj = 1 %[1 2 3 4] 
    subject = subjects{subj};
    
    try
        if contains(subject, 'E')
            folder = sprintf('../../out/%s', subject);
            filename = sprintf('Prediction_%s_%s-processed-rippA-0_Prewin-150_win.mat', subject, state);
            load(fullfile(folder, filename), 'nullStats')
            rippAstats = nullStats;
            
            filename = sprintf('Prediction_%s_%s-processed-0_Prewin-150_win.mat', subject, state);
            load(fullfile(folder, filename), 'nullStats', 'rippStats')

            
        elseif contains(subject, 'B1')
            if contains(subject,'dual')
                arrayConfig = 'lateral-medial';
                subject = 'B1';
            else
                arrayConfig = 'medial';
            end
            folder = sprintf('../../out/%s/%s', subject, arrayConfig);
            filename = sprintf('Prediction_%s_%s-processed-rippA-0_Prewin-150_win.mat', subject, state);
            load(fullfile(folder, filename), 'nullStats')
            rippAstats = nullStats;

            filename = sprintf('Prediction_%s_%s-processed-0_Prewin-150_win.mat', subject, state);
            load(fullfile(folder, filename))

            
        end
    catch 
        warning('could not load %s data for %s', state, subject)
        continue
    end

    coRippleRate = rippStats.coRippRate;
    channelDistances = rippStats.Distance;
    channelDistances(channelDistances == 0) = nan;
    unitDistances = nan(size(rippStats.predict));
    unitCoRip = nan(size(rippStats.predict));
    for iiA = 1:size(unitDistances,1)
            for iiB = 1:size(unitDistances,2)
                chA = rippStats.chA(iiA,iiB);
                chB = rippStats.chB(iiA,iiB);
                if chB > 96 
                    unitDistances(iiA, iiB) = 10;
                    unitCoRip(iiA, iiB) = coRippleRate(chA,chB-96);

                elseif chA > 0 && chB > 0
                    unitDistances(iiA, iiB) = channelDistances(chA,chB);
                    unitCoRip(iiA, iiB) = coRippleRate(chA,chB);
               
                end
            end
    end

    unitDistances = unitDistances(:);
    unitCoRip = unitCoRip(:);


    Ar = cellfun(@(X) isempty(X), rippStats.predict);
    Ara = cellfun(@(X) isempty(X), rippAstats.predict);
    An = cellfun(@(X) isempty(X), nullStats.predict);

    rippPredictAll = cell2mat(rippStats.predict(~Ar | ~An | ~Ara));
    rippAPredictAll = cell2mat(rippAstats.predict(~Ar | ~An | ~Ara));
    nullPredictAll = cell2mat(nullStats.predict(~Ar | ~An | ~Ara));

    unitDistances(Ar(:) | An(:) | Ara(:)) = [];
    unitCoRip(Ar(:) | An(:) | Ara(:)) = [];
    art = rippPredictAll > 3;

    rippPredictAll(art) = [];
    rippAPredictAll(art) = [];
    unitDistances(art) = [];
    nullPredictAll(art) = [];
    unitCoRip(art) = [];
    
    bins = quantile(unitDistances,8); %:0.01:0.9;
    bins(2:end+1) = bins;
    bins(1) = 0;
    bins(end+1) = 11;

    rippPredictAllBinned = [];
    rippAPredictAllBinned = [];
    nullPredictAllBinned = [];
    rippPredictAllBinnedSTD = [];
    rippAPredictAllBinnedSTD = [];
    nullPredictAllBinnedSTD = [];  
    coRippleBinned = [];
    coRippleBinnedSEM = [];
    
        
    clc
    nZeroPairs = [];
    c = 1;
    h =[];
    for iB = 1:(length(bins)-1)
        ii = unitDistances >= bins(iB) & unitDistances < bins(iB+1);
        sum(ii)
        rippPredictAllBinned = [rippPredictAllBinned, mean(rippPredictAll(ii))];
        rippAPredictAllBinned = [rippAPredictAllBinned, mean(rippAPredictAll(ii))];
        nullPredictAllBinned = [nullPredictAllBinned, mean(nullPredictAll(ii))];
        
        rippPredictAllBinnedSTD = [rippPredictAllBinnedSTD, std(rippPredictAll(ii))/ sqrt(sum(ii))];
        rippAPredictAllBinnedSTD = [rippAPredictAllBinnedSTD, std(rippAPredictAll(ii))/ sqrt(sum(ii))];
        nullPredictAllBinnedSTD = [nullPredictAllBinnedSTD, std(nullPredictAll(ii))/ sqrt(sum(ii))];
        

        coRippleBinned = [coRippleBinned, mean(unitCoRip(ii))];
        coRippleBinnedSEM = [coRippleBinnedSEM, std(unitCoRip(ii))/sqrt(sum(ii))];
    
        nZero = sum(rippPredictAll(ii) == 0 & nullPredictAll(ii) == 0) / length(ii);
        nZeroPairs = [nZeroPairs, nZero];
        
        [h(c),p,ci] = ttest(rippPredictAll(ii), nullPredictAll(ii));
                    
        c = c+1;
        
    end
    
    
    if subj == 2
        figure;
        hDistance = gcf;
        hDistance.Position = [990 458 126 253];
       
        nullPredictAll = cell2mat(nullStats.predict(:));
        pl1 = errorbar(1, mean(nullPredictAll), ...
                                 std(nullPredictAll)/sqrt(length(nullPredictAll)), ...
                                 'o'); hold on;
        pl1.Color = 0.8*clr(subj-1,:);
        pl1.LineWidth = 1.0;
        pl1.MarkerFaceColor = clr2(1,:);

        rippAPredictAll = cell2mat(rippAstats.predict(:));
        pl1 = errorbar(1, mean(rippAPredictAll), ...
                                 std(rippAPredictAll)/sqrt(length(rippAPredictAll)), ...
                                 'o'); hold on;
        pl1.Color = 0.8*clr(subj-1,:);
        pl1.LineWidth = 1.0;
        pl1.MarkerFaceColor = 'k'; clr2(1,:);

        rippPredictAll = cell2mat(rippStats.predict(:));
        pl2 = errorbar(1, mean(rippPredictAll), ...
                                 std(rippPredictAll)/sqrt(length(rippPredictAll)), ...
                                 'o'); hold on;
        pl2.Color = clr(subj-1,:);
        pl2.LineWidth = 1.0;
        pl2.MarkerFaceColor = clr2(8,:);

        
        ax = gca;
        maxLim = medialLim(2);
        minLim = -0.02;
        ylim([minLim maxLim ])
        ax.YTick = minLim:0.05:maxLim;
        ax.YTick = linspace(0, maxLim, 4);

       
        fig = gcf;
        fig.Color = 'w';
        box off 
        fname = sprintf('%s_%s_PredictDistance_SmoothWindow_%i_PSTHwindow_%i.pdf',[subject, 'dual'],state, smoothingWindow, win);
        savepdf(gcf, fullfile(exportDirec,fname));
    else
    
        hDistance = figure;
%         
        hDistance.Position = [995 463 230 253];
%         
        xVal = movmean(bins, 2) * 400;
        xVal(1) = [];

        pl1 = errorbar(xVal, nullPredictAllBinned, nullPredictAllBinnedSTD, 'o'); hold on;
        pl1.Color = 0.8*clr(subj,:);
        pl1.LineWidth = 1.0;
        pl1.MarkerFaceColor = clr2(1,:);

        idx = ~isnan(nullPredictAllBinned);
        p = polyfit(xVal(idx), nullPredictAllBinned(idx),1); %linear fit
        x1 = linspace(min(xVal),max(xVal));
        y1 = polyval(p,x1);
        l1 = plot(x1,y1); hold on;
        l1.Color =  clr2(1,:);

        pl1 = errorbar(xVal, rippAPredictAllBinned, rippAPredictAllBinnedSTD, 'o'); hold on;
        pl1.Color = 0.8*clr(subj,:);
        pl1.LineWidth = 1.0;
        pl1.MarkerFaceColor = 'k'; 

        idx = ~isnan(rippAPredictAllBinned);
        p = polyfit(xVal(idx), rippAPredictAllBinned(idx),1); %linear fit
        x1 = linspace(min(xVal),max(xVal));
        y1 = polyval(p,x1);
        l1 = plot(x1,y1); hold on;
        l1.Color =  'k'; clr2(1,:);
        
        
        pl2 = errorbar(xVal, rippPredictAllBinned,rippPredictAllBinnedSTD, 'o'); hold on;
        pl2.Color = clr(subj,:);
        pl2.LineWidth = 1.0;
        pl2.MarkerFaceColor = clr2(8,:);
        idx = ~isnan(rippPredictAllBinned);
        p = polyfit(xVal(idx), rippPredictAllBinned(idx),1); %linear fit
        x1 = linspace(min(xVal),max(xVal));
        y1 = polyval(p,x1);
        l1 = plot(x1,y1); hold on;
        l1.Color =  clr2(8,:);
        
        ax = gca;
        fig = gcf;
        fig.Color = 'w';
        maxLim = ax.YLim(2); 
        minLim = ax.YLim(1); 
        ylim([minLim maxLim ])
        ax.YTick = minLim:0.05:maxLim;

        if subj == 1
            medialLim = [minLim maxLim];
            ylim([-0.02 0.12 ])
            ax.YTick = linspace(0, 0.12, 4);

        end
        xlim([0 4500])

        box off

        fname = sprintf('%s_%s_PredictDistance_SmoothWindow_%i_PSTHwindow_%i.pdf',subject,state, smoothingWindow, win);
        savepdf(gcf, fullfile(exportDirec,fname));
    end

end


%% Predictability vs Cell Type (Fig. 4egi)
interactions = {'int->int','pyr->int', ...
                'int->pyr','pyr->pyr'};
close all
subject = 'B1'; %E1
folder = sprintf('../../out/%s/lateral-medial', subject);
% folder = sprintf('../../out/%s/', subject);
filename = sprintf('Prediction_%s_%s-processed-0_Prewin-150_win.mat', subject, state);
load(fullfile(folder, filename))
nullFiringRate = nullStats.FiringRate(:);
nullFiringRate(nullFiringRate == 0) = [];

rippFiringRate = rippStats.FiringRate(:);
rippFiringRate(rippFiringRate == 0) = [];
coRippleRate = rippStats.coRippRate;
channelDistances = rippStats.Distance;
channelDistances(channelDistances == 0) = nan;
unitDistances = nan(size(rippStats.predict));
unitCoRip = nan(size(rippStats.predict));
for iiA = 1:size(unitDistances,1)
        for iiB = 1:size(unitDistances,2)
            chA = rippStats.chA(iiA,iiB);
            chB = rippStats.chB(iiA,iiB);
            if chB > 96 
                unitDistances(iiA, iiB) = 10;
                unitCoRip(iiA, iiB) = coRippleRate(chA,chB-96);

            elseif chA > 0 && chB > 0
                unitDistances(iiA, iiB) = channelDistances(chA,chB);
                unitCoRip(iiA, iiB) = coRippleRate(chA,chB);
           
            end
        end
end

unitDistances = unitDistances(:);
unitCoRip = unitCoRip(:);


Ar = cellfun(@(X) isempty(X), rippStats.predict);
An = cellfun(@(X) isempty(X), nullStats.predict);

unitDistances(Ar(:) | An(:)) = [];

rippPredictAll = cell2mat(rippStats.predict(:));
nullPredictAll = cell2mat(nullStats.predict(:));

interactionType = rippStats.CellType(:);
interactionType(Ar(:) | An(:)) = [];

for ii = 1:length(interactions)

    iType = interactions{ii};
    iiType = cellfun(@(x) strcmp(x,iType), interactionType, 'UniformOutput', false);
    iiType = cell2mat(iiType);
    
    rippPredictType = rippPredictAll(iiType);
    nullPredictType = nullPredictAll(iiType);
    distancesType   = unitDistances(iiType);

    xAll = [ ]; 
    yAll = [ ];    
    xAllR = [ ]; 
    yAllR = [ ];    
    frAllR = [ ];    
    xAllN = [ ]; 
    yAllN = [ ];
    frAllN = [ ];    

    distAll = [];

    h2 = figure('Position', [1290 734 400 400]);
    for p = 1:length(rippPredictType)
        y = rippPredictType(p);
        x = nullPredictType(p);
        yFR = rippFiringRate(p);
        xFR = nullFiringRate(p);

        

        if (y - x) > 0  %&& y > 0 && x > 0
            f = distancesType(p)/max(unitDistances);  
            l = abs(y-x);
            c = 1 / (1 + exp(-l));

            xAll = [xAll x];
            yAll = [yAll y];
            xAllR = [xAllR x];
            yAllR = [yAllR y];
            frAllR = [ frAllR yFR];
        elseif (y - x) < -0  %&& y > 0 && x > 0
            f = distancesType(p)/max(unitDistances);
            l = abs(y-x);
            c = 1 / (1 + exp(-l));
            

            xAll = [xAll x];
            yAll = [yAll y];

            distAll = [distAll distancesType(p)];
            xAllN = [xAllN x];
            yAllN = [yAllN y];

            frAllN = [ frAllN xFR];

        else
%             pl = plot(x, y, '.', 'Color', [0.4 0.4 0.4]); hold on;
%             pl.MarkerSize = 5;
        end
    end
    dscatter(xAll', yAll'); hold on;                

    colormap(magma)
   

    pl = plot([-1.0 11], [-1.0 11], 'r--'); hold on;
    pl.LineWidth = 1.5;
    meanR = mean(yAll)
    semR = std(yAll); 
    meanN = mean(xAll);
    semN = std(xAll); 

    pl = plot([meanN meanN], [meanR-semR meanR+semR],'k-'); hold on;
    pl.LineWidth = 1.5;
    pl.Color = [0.5 0.5 0.5];
    pl = plot([meanN+semN meanN-semN], [meanR meanR],'k-'); hold on;
    pl.Color = [0.5 0.5 0.5];
    pl.LineWidth = 1.5;

    pl = plot(meanN,meanR,'dk');
    pl.MarkerFaceColor = [0.5 0.5 0.5]  ;  
    pl.MarkerSize = 8;
    ylabel('Co-Ripple Predictability')
    xlabel('No-Ripple Predictability')


    
    fprintf('%s   ', iType)

    [h p ci] = ttest(yAll, xAll, 'tail', 'right');
    fprintf('sigdiff h: %0.3f     p: %0.9f   \n\n', h, p);
    xlim([-0.5 1.0])
    ylim([-0.5 1.0])
    h2.Color = 'w';
    
    title(iType)

    fname = sprintf('%s_%s_PredictabilityScatter_SpikesOmitted_SmoothWindow_%i_PSTHwindow_%i_%s.pdf',subject,state,smoothingWindow, win,iType);
    savepdf(h2, fullfile(exportDirec,fname))
    
    figure('Position', [429 298 618 150]);    

    fprintf('%s   :', iType)

    
    [hst, bins] = histcounts(yAll - xAll, -1:0.1:1, 'Normalization', 'probability'); 
    xVal = movmean(bins, 2);
    xVal(1) = [];

    Y = [0 hst 0];
    X = [1 xVal xVal(end)];
    f = fill(X,Y, clr(2,:)); hold on;
    f.EdgeAlpha = 0;
    f.FaceAlpha = 0.5;


    pl = plot(xVal, hst); hold on;
    pl.LineWidth = 1.5;
    pl.Color = [207 178 132]/255;
    
    vline(0)
    xlim([-1.0 1.0])
    ylim([0 0.25])
    fig = gcf;
    fig.Color = 'w';
    title(iType)

    fname = sprintf('%s_%s_PredictabilityHist_SpikesOmitted_SmoothWindow_%i_PSTHwindow_%i_%s.pdf',subject, state,smoothingWindow, win,iType);
    savepdf(gcf, fullfile(exportDirec,fname));


end
%% %% Predictability vs Cell Type (Fig. 4fh)

close all
h =  findobj('type','figure');
nFig = length(h);
for subj = [1,3 4]
    subject = subjects{subj};
    
    fprintf('%s\n', subject);
    try
        if contains(subject, 'E')
            folder = sprintf('../../out/%s', subject);
            filename = sprintf('Prediction_%s_%s-processed-0_Prewin-150_win.mat', subject, state);
            load(fullfile(folder, filename))
        elseif contains(subject, 'B1')
            if contains(subject,'dual')
                arrayConfig = 'lateral-medial';
                subject = 'B1';
            else
                arrayConfig = 'medial';
            end
            folder = sprintf('../../out/%s/%s', subject, arrayConfig);
            filename = sprintf('Prediction_%s_%s-processed-0_Prewin-150_win.mat', subject, state);
            load(fullfile(folder, filename))
        end
    catch 

        warning('could not load %s data for %s', state, subject)
        continue
    end
    unitDistances = rippStats.Distance(:);

    nullFiringRate = nullStats.FiringRate(:);
    
    rippFiringRate = rippStats.FiringRate(:);
    
    Ar = cellfun(@(X) isempty(X), rippStats.predict);
    An = cellfun(@(X) isempty(X), nullStats.predict);

    rippPredictAll = cell2mat(rippStats.predict(:));
    nullPredictAll = cell2mat(nullStats.predict(:));

    interactionType = rippStats.CellType(:);
    interactionType(Ar | An) = [];
    predictTypeAll = cell(1,4);

    for ii = 1:length(interactions)
    
        iType = interactions{ii};
        iiType = cellfun(@(x) strcmp(x,iType), interactionType, 'UniformOutput', false);
        iiType = cell2mat(iiType);
        
        rippPredictType = rippPredictAll(iiType);
        nullPredictType = nullPredictAll(iiType);
    
        xAll = [ ]; 
        yAll = [ ];    
        xAllR = [ ]; 
        yAllR = [ ];    
        frAllR = [ ];    
        xAllN = [ ]; 
        yAllN = [ ];
        frAllN = [ ];    
    
        distAll = [];
    
        for p = 1:length(rippPredictType)
            y = rippPredictType(p);
            x = nullPredictType(p);
            yFR = rippFiringRate(p);
            xFR = nullFiringRate(p);
    
            
    
            if (y - x) > 0  %&& y > 0 && x > 0
                l = abs(y-x);
                c = 1 / (1 + exp(-l));
    
                xAll = [xAll x];
                yAll = [yAll y];
                xAllR = [xAllR x];
                yAllR = [yAllR y];
                frAllR = [ frAllR yFR];
                frAllN = [ frAllN xFR];
            elseif (y - x) < -0  %&& y > 0 && x > 0
                l = abs(y-x);
                c = 1 / (1 + exp(-l));
        
                
    
                xAll = [xAll x];
                yAll = [yAll y];
    
                xAllN = [xAllN x];
                yAllN = [yAllN y];
    
                frAllR = [ frAllR yFR];
                frAllN = [ frAllN xFR];    
            else
    %             pl = plot(x, y, '.', 'Color', [0.4 0.4 0.4]); hold on;
    %             pl.MarkerSize = 5;
            end
        end
        
        predictTypeAll{ii} = [predictTypeAll{ii} yAll];
        figure(nFig + ii); 
        set(gcf,'Position',[767 858 385 2*97])

        fprintf('%s   :', iType)
        
        [hst, bins] = histcounts(yAll - xAll, -2:0.1:2, 'Normalization', 'probability'); 
        xVal = movmean(bins, 2);
        xVal(1) = [];
    
        Y = [0 hst 0];
        X = [1 xVal xVal(end)];
        f = fill(X,Y, clr(subj,:)); hold on;
        f.EdgeAlpha = 0;
        f.FaceAlpha = 0.5;
    
    
        pl = plot(xVal, hst); hold on;
        pl.LineWidth = 1.5;
        pl.Color = clr(subj,:);
        vline(0)        
        xlim([-1.0 1.0])
        fig = gcf;
        fig.Color = 'w';
        box off                         
        [H,P,CI] = ttest(yAll - xAll);
        fprintf(' ttest p: %0.4f\n', P)
        title(iType)
        iType(iType == '>') = [];
        
        if subj == length(subjects)
            fname = sprintf('%s_PredictabilityHist_SpikesOmitted_SmoothWindow_%i_PSTHwindow_%i_%s.pdf',state,smoothingWindow, win,iType);
            savepdf(gcf, fullfile(exportDirec,fname));
        end
    
    end
    
    [h, p, ci, stats] = ttest2(predictTypeAll{4},predictTypeAll{1}, 'tail', 'right');
    fprintf('pyrpyr vs intint %.5f p = %.5f\n', mean(predictTypeAll{4}) - mean(predictTypeAll{1}), p)
end


