
close all
clear
clc

addpath(genpath('../../code'))


%%
subjects = {'B1', 'B1-dual', 'E1', 'E2'};
state = 'NREM';
clr = [237/255 110/255 87/255; %B1
       255/255 66/255 161/255; %B1-dual
       108/255 227/255 208/255; %E1
       129/255 213/255 83/255]; %E2

exportDirec = '../../figures';
if ~isfolder(exportDirec); mkdir(exportDirec); end

computeCoRip = false;
%% Predictability vs Distance
close all
for subj = 3
    subject = subjects{subj};
    
    try
        if contains(subject, 'E')
            folder = sprintf('../../out/%s', subject);
            filename = sprintf('Prediction_%s_%s-processed-0_Prewin-150_win.mat', subject, state);
            load(fullfile(folder, filename))
            
            goodChannels = 1:96;
            badChannels = find(~ismember(1:96, goodChannels));

        elseif contains(subject, 'B1')
            if contains(subject,'dual')
                arrayConfig = 'lateral-medial';
                subject = 'B1';
                chA_all = 1:96;
                chB_all = 97:192;
            else
                arrayConfig = 'medial';
                chA_all = 1:96;
                chB_all = 1:96;
            end
            folder = sprintf('../../out/%s/%s', subject, arrayConfig);
            filename = sprintf('Prediction_%s_%s-processed-0_Prewin-150_win.mat', subject, state);
            load(fullfile(folder, filename))

           
            
        end

        
        goodChannels = 1:96; %unique(cell2mat(units(:,1)));
        badChannels = find(~ismember(1:96, goodChannels));
    catch 
        warning('could not load %s data for %s', state, subject)
        continue
    end
    

    
    coRippleRate = rippStats.coRippRate; 
    channelDistances = rippStats.Distance;

    
    coRippleRate(badChannels, :) = [];
    coRippleRate(:, badChannels) = [];
    channelDistances(badChannels, :) = [];
    channelDistances(:, badChannels) = [];
    

    channelDistances(coRippleRate == 0) = [];
    coRippleRate(coRippleRate == 0) = [];

    channelDistances = channelDistances(:);
    coRippleRate = coRippleRate(:);
    
    rippPredictAll = cell2mat(rippStats.predict(:));
    nullPredictAll = cell2mat(nullStats.predict(:));
    
    bins  = linspace(400,4500, 15); %+(8*400);
    bins = bins/400;
    dW = diff(bins);
    bins(end) = [];
    
    rippPredictAllBinned = [];
    nullPredictAllBinned = [];
    rippPredictAllBinnedSTD = [];
    nullPredictAllBinnedSTD = [];   
    
    coRippleBinned = [];
    coRippleBinnedSTD = [];
        
    clc
    nZeroPairs = [];
    c = 1;
    h = [];
    for W = bins
        ii = channelDistances > W & channelDistances < W + dW(c); % & MannUAll == 1;
        sum(ii)

        coRippleBinned = [coRippleBinned, mean(coRippleRate(ii))];
        coRippleBinnedSTD = [coRippleBinnedSTD, std(coRippleRate(ii))/sqrt(sum(ii))];
            
        fprintf('%0.0f %0.0f\n', 400*W, 400*(W+dW(c)))
        
        c = c+1;
        
    end
    
    if subj == 2
        
        goodChannels = 1:96;
        badChannelsM = find(~ismember(1:96, goodChannels));

        goodChannels = 1:96;
        badChannelsL = [1:96, find(~ismember(1:96, goodChannels)) + 96];

        coRippleRate(badChannelsM, :) = [];
        coRippleRate(:, badChannelsL) = [];

        figure(3);
        hDistance = gcf;
        hDistance.Position = [189 947 65 250];
        
        pl1 = errorbar(1.5, mean(coRippleRate(:)),std(coRippleRate(:))/sqrt(length(coRippleRate(:))), 'o'); hold on;
        pl1.Color = clr(subj,:);
        pl1.LineWidth = 1.0;
        
        ax = gca;
        ylim([0.1 0.25])
        xlim([1 2])
    
        fig = gcf;
        fig.Color = 'w';

        box off

        fname = sprintf('%s_%s_coRipDistance_SmoothWindow_%i_PSTHwindow_%i.pdf',[subject, 'dual'],state, smoothingWindow, win);
        savepdf(gcf, fullfile(exportDirec,fname));
    elseif subj == 1
        figure(1);
        hDistance = gcf;
        hDistance.Position = [189 947 266 250];
        
        pl1 = errorbar(400*bins+200, coRippleBinned, coRippleBinnedSTD, 'o'); hold on;
        pl1.Color = clr(subj,:);
        pl1.LineWidth = 1.0;
        x = 400*bins+200;
        y = coRippleBinned;
        f = @(beta, x) beta(1)*exp(-beta(2)*x) + beta(3);
    
        % Define the error function to minimize (sum of squares)
        err = @(beta) sum((f(beta, x) - y).^2);
        
        % Use fminsearch to minimize the error function
        beta0 = [1, 0, mean(y)]; % Initial guess for parameters
        beta_hat = fminsearch(err, beta0);

        pl = plot(x, f(beta_hat, x)); hold on;
        pl.Color = clr(subj, :);
        
        ax = gca;
        ylim([0.10 0.25])
        box off
    
        fig = gcf;
        fig.Color = 'w';
        fname = sprintf('%s_%s_coRipDistance_SmoothWindow_%i_PSTHwindow_%i.pdf',subject,state, smoothingWindow, win);
        savepdf(gcf, fullfile(exportDirec,fname));
    else
        

        figure(2);
        hDistance = gcf;
        hDistance.Position = [189 947 266 250];
        
        pl1 = errorbar(400*bins+200, coRippleBinned, coRippleBinnedSTD, 'o'); hold on;
        pl1.Color = clr(subj,:);
        pl1.LineWidth = 1.0;
        x = 1:length(coRippleBinned);
        y = coRippleBinned;
        f = @(beta, x) beta(1)*exp(-beta(2)*x) + beta(3);
    
        % Define the error function to minimize (sum of squares)
        err = @(beta) sum((f(beta, x) - y).^2);
        
        % Use fminsearch to minimize the error function
        beta0 = [1, 1, mean(y)]; % Initial guess for parameters
        beta_hat = fminsearch(err, beta0);

        pl = plot(400*bins+200, f(beta_hat, x)); hold on;
        pl.Color = clr(subj, :);
        ax = gca;
        if strcmp(state, 'NREM')
%             ylim([0.10 0.25])
        else
            ylim([0.0 0.15])
        end

        box off
    
        fig = gcf;
        fig.Color = 'w';
    end
    
    if subj == 3
        fname = sprintf('%s_%s_coRipDistance_SmoothWindow_%i_PSTHwindow_%i.pdf','All',state, smoothingWindow, win);
        savepdf(gcf, fullfile(exportDirec,fname));
    end

end




%%











