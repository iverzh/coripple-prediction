
close all 
clear 
clc


%% Upstate ripple co fire
subjects = {'B1', 'B1-dual', 'E1', 'E2'};
state = 'NREM';
clr = [237/255 110/255 87/255; %B1
       255/255 66/255 161/255; %B1-dual
       108/255 227/255 208/255; %E1
       129/255 213/255 83/255]; %E2

exportDirec = '../../figures';
figure('Position', [857 687 461 413]);

nullAll = [];
USUSAll = [];
ripprippAll = [];   
USrippAll = [];

for subj = 1:4  
    subject = subjects{subj};
    exportSO = '../../out/cofire';
    filename = sprintf('%s_SO_CoFire.mat', subject);
    load(fullfile(exportSO,filename))
    

    keep = nullCoFire(:) > 0 & USUSCoFire(:) > 0 & ripprippCoFire(:) > 0 & USrippCoFire(:) > 0;
    nullAll = [nullAll nullCoFire(keep)'];
    USUSAll = [USUSAll USUSCoFire(keep)'];
    ripprippAll = [ripprippAll ripprippCoFire(keep)'];
    USrippAll = [USrippAll USrippCoFire(keep)'];
    

    nullPlot = mean(nullCoFire(:), 'omitnan');
    USUSPlot = mean(USUSCoFire(:), 'omitnan') / nullPlot;
    ripprippPlot = mean(ripprippCoFire(:), 'omitnan') / nullPlot;
    USrippPlot = mean(USrippCoFire(:), 'omitnan') / nullPlot;
    
    USUSsem = computeSEM(USUSCoFire / nullPlot, 1);
    ripprippsem = computeSEM(ripprippCoFire / nullPlot, 1);
    USrippsem= computeSEM(USrippCoFire / nullPlot, 1);  
    
    er = errorbar(1:3, [USUSPlot, ripprippPlot, USrippPlot] * 100, [USUSsem, ripprippsem, USrippsem] * 100, 'o-'); hold on;
    er.Color = clr(subj,:);
    er.MarkerFaceColor = er.Color;
    er.LineWidth = 1.5;
end
xlim([0.5 3.5])
ylim([0 900])
hline(100); hold on;
ylabel('CoFire Event Rate [% of control]')
ax = gca;
preformatLabels = {'co-US', 'ripple (-)', ...
                   'US (-)', 'co-ripple', ...
                   'co-US', 'co-ripple'};
% formatStr = repmat('%s\\newline%s\n', [1 length(preformatLabels)/2])
tickLabels = strtrim(sprintf('%s\\newline%s\n', preformatLabels{:}));
ax.XTick = 1:3;
ax.XTickLabel = tickLabels;
% ax.I
% title([subject, ' across'])
box off
% l = legend(subjects(1:end), 'Location', 'northwest');
grid on
fig = gcf;
fig. Color = 'w';
title('Co-firing with 10 ms')
savepdf(gcf, fullfile(exportDirec, 'SlowWaveRippleCoFire_10ms.pdf'))


[h,p] = signrank(USUSAll, nullAll, 'tail', 'right');
fprintf('USUS vs null: %.4f\n', h)
[h,p] = signrank(ripprippAll, nullAll, 'tail', 'right');
fprintf('rippripp vs null: %.4f\n', h)
[h,p] = signrank(USrippAll, nullAll, 'tail', 'right');
fprintf('rippUS vs null: %.4f\n', h)

[h,p] = signrank(USrippAll, USUSAll, 'tail', 'right');
fprintf('rippUS vs USUS: %.4f\n', h)

[h,p] = signrank(USrippAll, ripprippAll, 'tail', 'right');
fprintf('rippUS vs rippripp: %.4f\n', h)

[h,p] = signrank(ripprippAll, USUSAll, 'tail', 'right');
fprintf('rippripp vs USUS: %.4f\n', h)

mu = mean(USrippAll, 'omitnan') * 1e3;
sem = computeSEM(USrippAll ,1) * 1e3;
fprintf('USrippAll: %.4f %s %.4f \n', mu, char(177), sem)

mu = mean(ripprippAll, 'omitnan') * 1e3;
sem = computeSEM(ripprippAll ,1) * 1e3;
fprintf('ripprippAll: %.4f %s %.4f \n', mu, char(177), sem)

mu = mean(USUSAll, 'omitnan') * 1e3;
sem = computeSEM(USUSAll ,1) * 1e3;
fprintf('USUSAll: %.4f %s %.4f \n', mu, char(177), sem)

mu = mean(nullAll, 'omitnan') * 1e3;
sem = computeSEM(nullAll ,1) * 1e3;
fprintf('nullAll: %.4f %s %.4f \n', mu, char(177), sem)

%%
close all
direc =   '../../out/cofire';
state = 'wake';
figure('Position', [852 682 343 406]);
nCo = [];
nNo = [];
durationAll = [];
for subj = [1 3] %[1 3 4] [1 3]

    filename = sprintf('%s_coFireCoRipplePeriod_%s.mat', subjects{subj}, state);
    load(fullfile(direc,filename))
   
   
    for uA = 1:size(coRipDur, 1)
        for uB = 1:size(coRipDur, 2)
            coRip = coRipSpike{uA,uB};
            noRip = noRipSpike{uA,uB};
            dur = sum(coRipDur{uA,uB}) / 1e3;

            if (~isempty(coRip) || ~isempty(noRip)) %&& strcmp(interactionType{uA, uB}, interact) % || strcmp(interactionType{uA, uB}, 'int-pyr'))
                sp = coRip(:,coRip(1,:) > 0 & coRip(2,:) > 0);
                nCoSpikeCo = size(sp,2); %/dur; % / (dur);
                
                sp = noRip(:,noRip(1,:) > 0 & noRip(2,:) > 0);
                nCoSpikeNo = size(sp,2); %/dur; % / (dur);
%                     dur = dur + (sum(coRipDur{uA,uB}) / 1e3);
                
                nCo = [nCo, nCoSpikeCo];
                nNo = [nNo, nCoSpikeNo];
                durationAll = [durationAll, dur];
            end
        end
    end
    
    if ~contains(subjects{subj}, 'dual')
        nCoSpikeCoDist = [];
        nCoSpikeNoDist = [];
        nCoSpikeCoDistSEM = [];
        nCoSpikeNoDistSEM = [];
        bins = quantile(distanceAll, 0:1/7:1); %:0.01:0.9;
        
        zAll = zscore([nCoSpikeCoAll nCoSpikeNoAll]);
        nCoSpikeCoAll = zAll(1:length(nCoSpikeCoAll));
        nCoSpikeNoAll = zAll(length(nCoSpikeCoAll)+1:end);

        for iB = 1:(length(bins)-1)
            ii = distanceAll >= bins(iB) & distanceAll < bins(iB+1);
            
            nCoSpikeCoDist = [nCoSpikeCoDist, mean(nCoSpikeCoAll(ii))];
            nCoSpikeNoDist = [nCoSpikeNoDist, mean(nCoSpikeNoAll(ii))];

            nCoSpikeCoDistSEM = [nCoSpikeCoDistSEM, computeSEM(nCoSpikeCoAll(ii), 1)];
            nCoSpikeNoDistSEM = [nCoSpikeNoDistSEM, computeSEM(nCoSpikeNoAll(ii),1)];
        
        
        end

        xVal = movmean(bins, 2) * 400;
        xVal(1) = [];

        
        
        clr2 = brewermap(10,'Paired');
        
        p = polyfit(xVal, nCoSpikeCoDist,1); %linear fit
        x1 = linspace(min(xVal),max(xVal));
        y1 = polyval(p,x1);
%         l1 = plot(x1,y1); hold on;
%         l1.Color =  clr2(8,:);
%         l1.Color =  clr(subj,:);
        b1 = errorbar(xVal, nCoSpikeCoDist, nCoSpikeCoDistSEM, 'o-'); hold on;
        b1.MarkerFaceColor = clr(subj,:);
%         b1.Color = clr2(8,:);
        b1.Color = clr(subj,:);
        b1.LineWidth = 1.0;
        
        
        
        
        p = polyfit(xVal, nCoSpikeNoDist,1); %linear fit
        x1 = linspace(min(xVal),max(xVal));
        y1 = polyval(p,x1);
%         l2 = plot(x1,y1); hold on;
%         l2.Color =  clr2(1,:);
        l2.Color =  0.6*clr(subj,:);
        
        b2 = errorbar(xVal, nCoSpikeNoDist, nCoSpikeNoDistSEM, 'o--'); hold on;
        b2.MarkerFaceColor = 0.6*clr(subj,:);
        b2.Color = 0.6*clr(subj,:);
        b2.LineWidth = 1.0;
        b2.MarkerSize = 5;
        
        
        ylabel('CoFire Event Rate [z-score]')
        % legend('CoRipples','Control')
        
        % yyaxis right
        % pl = plot(nCoSpikeCoDist./nCoSpikeNoDist);
        % pl.Color = [180/255 160/255 101/255];
        % pl.LineWidth = 2;
        % title(state)
        fig = gcf;
        fig.Color = 'w';
        box off
        
        xlim([0 4000])
    
        ax = gca;
        ax.XTick = 0:500:4000;
%         ax.XTickLabel = ax.XTick * 500;
%         ax.YTick = linspace(-ax.YTick(end),ax.YTick(end),6);
%         ax.YTickLabel = arrayfun(@(X) sprintf('%1.1f', X), ax.YTick, 'UniformOutput', false);
    
%         savepdf(gcf, fullfile(exportDirec, [subjects{subj},'CoRippleCoFire_wake.pdf']))
    else
       
        figure('Position', [852.0000 681 149.0000 413]);
        
        zAll = zscore([nCoSpikeCoAll nCoSpikeNoAll]);
        nCoSpikeCoAll = zAll(1:length(nCoSpikeCoAll));
        nCoSpikeNoAll = zAll(length(nCoSpikeCoAll)+1:end);
        
        clr2 = brewermap(10,'Paired');
        
        
%         b1 = errorbar(1:length(nCoSpikeCoDist), nCoSpikeCoDist, nCoSpikeCoDistSEM, 'o'); hold on;
        b1 = plot(3, mean(nCoSpikeCoAll), 'o'); hold on;
        b1.MarkerFaceColor = clr(subj,:);
        b1.Color = clr(subj,:);
        b1.LineWidth = 1.0;
        
        
        
        
        
%         b2 = errorbar(1:length(nCoSpikeNoDist), nCoSpikeNoDist, nCoSpikeNoDistSEM, 'o'); hold on;
        b2 = plot(3, mean(nCoSpikeNoAll), 'o'); hold on;
        b2.MarkerFaceColor = 0.6*clr(subj,:);
        b2.Color = 0.6*clr(subj,:);
        b2.LineWidth = 1.0;
        
        
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
        
        xlim([2.5 3.5])
    
        

%         close
    end
    
    


end
% ylim([-0.4000 0.6]) %ax.YTick(end)])
ylim([-0.4000 1]) %ax.YTick(end)])

ax = gca;
%         ax.XTick = 0:2.5:10;
%         ax.XTickLabel = ax.XTick * 500;
% ax.YTick = linspace(ax.YTick(1),ax.YTick(end),14);
% ax.YTickLabel = arrayfun(@(X) sprintf('%1.1f', X), ax.YTick, 'UniformOutput', false);

savepdf(gcf, fullfile(exportDirec, [subjects{subj},'CoRippleCoFire_NREM.pdf']))
keep = nNo > 0 & nCo > 0; 
FRco = nCo(keep) ./ durationAll(keep);
fprintf('co ripple FR: %0.3f %s %0.3f ', mean(FRco), char(177), computeSEM(FRco,1))
FRno = nNo(keep) ./ durationAll(keep);
fprintf('no ripple FR: %0.3f %s %0.3f\n', mean(FRno), char(177), computeSEM(FRno,1))
% mean(FRco) / mean(FRno)
% [h, p, ci, stats] = signrank(FRco, FRno, 'tail', 'both');
% [h, p] = signrank(FRco, FRno, 'tail', 'both')

%%
close all
subjects = {'B1', 'B1-dual', 'E1', 'E2'};
state = 'wake';
clr = [237/255 110/255 87/255; %B1
       255/255 66/255 161/255; %B1-dual
       108/255 227/255 208/255; %E1
       129/255 213/255 83/255]; %E2
ShufflePValues = cell(1,4);
for subj = [1 3] %[2 3 4]
    subject = subjects{subj};

    load(sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/CoFire/%s_%s_coFire-v2.mat', subject, state))
    load(sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/CoFire/%s_%s_coFire.mat', subject, state), 'interactionType')
%         load(sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/CoFire/%s_%s_coFire_allms.mat', subject, state))
    load(sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/CoFire/%s_%s_coFire_noRip_allms.mat', subject, state))
    types = {'pyr-pyr', 'pyr-int', 'int-pyr', 'int-int'};
    for ii = 1:length(types)
       
%         load(sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/CoFire/%s_coFireCoRipplePeriod_%s.mat',subject, state), 'coRipDur')
        % boxplot(coFireProbCoR_2(coFireProbCoR_2 > 0), 'BoxStyle','filled');
    %     subplot(2,2,ii)
        ind = strcmp(interactionType,types(ii));
%         checkEmpty = cellfun(@(X) isempty(X), coFireProbCoR);
        for iiA = 1:size(coFireProbCoR,1)
            for iiB = 1:size(coFireProbCoR,2)
                if isempty(coFireProbCoR{iiA,iiB})
                    coFireProbCoR{iiA,iiB} = zeros(1,101);
                end
            end
        end

         for iiA = 1:size(coFireProbNoR,1)
            for iiB = 1:size(coFireProbNoR,2)
                if isempty(coFireProbNoR{iiA,iiB})
                    coFireProbNoR{iiA,iiB} = zeros(1,101);
                end
            end
        end
    %     ind = ind & coFireProbCoR > 10;
        ccg = cell2mat(coFireProbCoR(ind));
        ccgNoR = cell2mat(coFireProbNoR(ind));

%         coRippleDur = cellfun(@(X) sum(X), coRipDur);
        coRippleDurInteract = coRippleDur(ind); 
%         FRs = coRippleFRa(ind) < 2 & coRippleFRb < 2;
        ind2 = coRippleFRa(ind) < inf & coRippleFRb(ind) < inf & coRippleFRa(ind) > 1 & coRippleFRb(ind) > 1 & (sum(ccg,2) > 50); %& FRs < 1;
%         ind2 = sum(ccg, 2) > 0;
        ccg = ccg(ind2,:);
        sum(ind2)
       
        ccgNoR = ccgNoR(ind2,:);    %     ccg = smoothdata(ccg, 2, 'gaussian', 5);
        ccg_shuff_cell = coFireProbCoR_Shuff(ind);
        ccg_shuff_cell = ccg_shuff_cell(ind2,:);
        ccg_shuff = nan(length(ccg_shuff_cell), 1000, 101);
        for iS = 1:length(ccg_shuff_cell)
            ccg_shuff(iS,:,:) = ccg_shuff_cell{iS};
        end
        
        coF = 10;
        N =[];
        for iP = 1:size(ccg_shuff,1)
            ccgShuffPair = squeeze(ccg_shuff(iP,:,:)); 
            coShuff = sum(ccgShuffPair(:,51-coF:51+coF), 2);
            coReal = sum(ccg(iP,51-coF:51+coF));
            coRealNoR = sum(ccgNoR(iP,51-coF:51+coF));
            N(iP) = sum(coShuff < coReal);
%     
            figure(9);
            subplot(2,1,1)
            histogram(coShuff); hold on; vline(coReal); hold on; vline(coRealNoR)
            subplot(2,1,2)
            plot(-50:50, ccg(iP,:)); hold on;
            plot(-50:50, ccgNoR(iP,:)); hold on;
            titleStr = sprintf('%i / %i  p: %.4f', sum(ccg(iP, 26:50)), sum(ccg(iP, [26:50,51:76])), sum(coShuff < coReal)/1000);
            title(titleStr)
            waitforbuttonpress; clf;
%      
    
        
        end
    
        N = 1 - N/1000;

        ShufflePValues{ii} = [ShufflePValues{ii} N];
%         ccg_shuff = cell2mat(coFireProbCoR_Shuff(ind));
%         ccg_shuff = reshape(ccg_shuff, [], 100, 101);
%         ccg_shuff = ccg_shuff(ind2,:,:);
        

%         ccg_shuff_mu = squeeze(mean(ccg_shuff, 2));
%     %     ccg_shuff_mu = smoothdata(ccg_shuff_mu, 2, 'gaussian', 5);
% 
%         coRippleDur = coRippleDur(ind2);
%     
%         coFire = nan(size(ccg,1), 50);
%         for coF = 1:50; coFire(:,coF) = sum(ccg(:,51-coF:51+coF),2); end
% 
%         coFireNoR = nan(size(ccg,1), 50);
%         for coF = 1:50; coFireNoR(:,coF) = sum(ccgNoR(:,51-coF:51+coF),2); end
%         
%         coFireShuff = nan(size(ccg_shuff,1), size(ccg_shuff,2), 50);
%         for coF = 1:50; for s = 1:100; coFireShuff(:,s,coF) = sum(squeeze(ccg_shuff(:,s,51-coF:51+coF)),2); end; end
% %         coFireAna = cell2mat(coFireProbCoR_Ana(ind));
% %         coFireAna = coFireAna(ind2, :);
%         
%         
%         
%         coFireShuff_mu = squeeze(mean(coFireShuff, 2));
%         
%         % figure; 
%         % yyaxis left
%         % % plot(mean(coFire)./mean(coFireShuff_mu)); hold on;
%         % plot(mean(coFire), 'r-'); hold on;
%         % plot(mean(coFireShuff_mu), 'k-'); hold on;
%         pAll = [];
%         for i = 1:50
%             [h,p] = ttest(coFire(:,i) ./ coRippleDur, coFireShuff_mu(:,i) ./ coRippleDur, 'tail', 'right');
%             pAll(i) = p;
%         end
%     
%         pMask = pAll < 0.05;
%         b = mask2bounds(pMask);
%         
%         fprintf('%s: %i ms\n', subject, sum(pMask))
%     
%     
%     
%     %     bar([mean( (coFireProbCoR(ind) - coFireProbCoR_Shuff(ind))/mean(coFireProbCoR_Shuff(ind)))]); hold on;
%     
%         figure('Position',[872 823 278 301]);
% %         mu = mean((coFire - coFireShuff_mu)./coFireShuff_mu); 
% %         sigma =  std(coFire./coFireShuff_mu) / sqrt(size(coFire,1));
%         mu = mean(coFire ./ coRippleDur) * 1000;
%         sigma =  std(coFire ./ coRippleDur) / sqrt(size(coFire,1));
%         [bl , ba] = boundedline(1:50, mu, sigma, '-'); hold on;
%         ba.FaceColor = clr(subj,:);
%         ba.FaceAlpha = 0.5;
%         bl.Color = clr(subj,:);
%         bl.LineWidth = 2;
% 
%         mu = mean(coFireShuff_mu ./ coRippleDur)* 1000;
%         sigma =  std(coFireShuff_mu ./ coRippleDur) / sqrt(size(coFireShuff_mu,1));
% %         mu = mean(coFireAna);
%         [bl , ba] = boundedline(1:50, mu, sigma, '-'); hold on;
%         ba.FaceColor = 0.8*clr(subj,:);
%         ba.FaceAlpha = 0.5;
%         bl.Color = 0.8*clr(subj,:);
%         bl.LineWidth = 2;
% 
%         mu = mean(coFireNoR ./ coRippleDur) * 1000;
%         sigma =  std(coFireNoR ./ coRippleDur) / sqrt(size(coFireNoR,1));
%         [bl , ba] = boundedline(1:50, mu, sigma, '-'); hold on;
%         ba.FaceColor = 0.6*clr(subj,:);
%         ba.FaceAlpha = 0.5;
%         bl.Color = 0.6*clr(subj,:);
%         bl.LineWidth = 2;
% 
% %         sigTimes = b(1):b(end);
% %         pl = plot(sigTimes, repmat(max(mu+sigma), [1 length(sigTimes)]));
% %         pl.Color = [207 181 59]/255;
% %         pl.LineWidth = 3;
%     
%     %     histogram(coFireProbCoR(ind) - coFireProbCoR_Shuff(ind), 100)
%             
% 
% 
%        
%       
%     
%     %     histogram(coFireProbCoR(ind) - coFireProbCoR_Shuff(ind), 100)
%             
%         h = hline(0);
%         h.Color = 'k';
%         h.LineStyle = '-';
%     
%         box off;
%         grid on
%     
%         fig = gcf;
%         fig.Color = 'w';
%         xlim([0 50])
%         ax = gca;
%         ax.XTick = 0:10:50;
%     %     ylim([0 3])  
% %         savepdf(gcf, fullfile(exportDirec, sprintf('%s_coFireShuff_wake.pdf', subject)))
% 
%         figure('Position',[435 638 178 143]);
%         pl = plot(-50:50, mean(ccg)); hold on; 
%         pl.Color = clr(subj,:);
%         pl.LineWidth = 2.5;
%         pl = plot(-50:50, mean(ccgNoR)); hold on; 
%         pl.Color = 0.6*clr(subj,:);
%         pl.LineWidth = 2.5;
%         pl = plot(-50:50, mean(ccg_shuff_mu)); hold on;
%         pl.Color = 0.8*clr(subj,:);
%         pl.LineWidth = 2.5;
%         vline([-33:11:33])
% 
%         fig = gcf;
%         fig.Color = 'w';
% %         savepdf(gcf, fullfile(exportDirec, sprintf('%s_CCG_wake.pdf', subject)))
% 
%         ccgSmooth = smoothdata(ccg, 2, 'gaussian', 10);
%         win = 45; %crop window to ignore edge effects
%         winSamples = [-win:win] + 51; 
%         ccgPeaks = nan(1,length(ccg));
%         ccgBinom = nan(1,length(ccg));
%         for iC = 1:length(ccgPeaks)
%             peak = find(ccgSmooth(iC,winSamples) == max(ccgSmooth(iC,winSamples))) - round(length(winSamples)/2);
%             if length(peak) < 2
%                 ccgPeaks(iC) = peak(abs(peak) == max(abs(peak))); %keep closest CCG value to 0
%             end
% 
%             s = sum(ccg(iC, 26:50));
%             n = sum(ccg(iC, [26:50, 52:76]));
%             ccgBinom(iC) = myBinomTest(s,n,0.5);
%         end
% 
%         figure; 
%         histogram(ccgPeaks, -50:50);



    end

%     savepdf(gcf, fullfile(exportDirec, sprintf('%s_coFireShuff.pdf', subject)))
end   

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(ShufflePValues{1},.05,'pdep','yes');
sum(h) / length(h)
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh([ShufflePValues{2} ShufflePValues{3}] ,.05,'pdep','yes');
sum(h) / length(h)

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(ShufflePValues{4},.05,'pdep','yes');
sum(h) / length(h)

%% 
figure(9)
figure(10);
pVals = [];
for coF = 10 %1:50
    N =[];
    for iP = 1:length(ccg_shuff)
        ccgShuffPair = squeeze(ccg_shuff(iP,:,:)); 
        coShuff = sum(ccgShuffPair(:,51-coF:51+coF), 2);
        coReal = sum(ccg(iP,51-coF:51+coF));
        N(iP) = sum(coShuff < coReal);

        figure(9);
        subplot(2,1,1)
        histogram(coShuff); hold on; vline(coReal);
        subplot(2,1,2)
        plot(-50:50, smoothdata(ccg(iP,:), 'gaussian', 5));
        titleStr = sprintf('%i / %i  p: %.4f', sum(ccg(iP, 26:50)), sum(ccg(iP, [26:50,51:76])), ccgBinom(iP));
        title(titleStr)
        waitforbuttonpress; clf;


    
    end
    figure(10);

    N = 1 - N/1000;
    pVals = [pVals; N];
    F = sum(N < 0.05)/length(ccg_shuff);
    plot(coF, F, 'o'); hold on;
end

% save([subject,'_NREM_CCG_Analysis.mat'], 'pVals')
%% co fire avalanche plot 
subjects = {'B1', 'B1-dual', 'E1', 'E2'};
state = 'wake';
clr = [237/255 110/255 87/255; %B1
       255/255 66/255 161/255; %B1-dual
       108/255 227/255 208/255; %E1
       129/255 213/255 83/255]; %E2
CoupledAllSubj = [];
unCoupledAllSubj = [];
nAllSubj = [];
for subj = [1 3]
    subject = subjects{subj};
    load(sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/CoFire/%s_%s_CoFireAvalanche.mat', subject, state))
    keep = ~isnan(CoupledAll) & isfinite(CoupledAll);
    CoupledAllSubj = [CoupledAllSubj CoupledAll(keep)];
    unCoupledAllSubj = [unCoupledAllSubj unCoupledAll(keep)];
    nAllSubj = [nAllSubj nAll(keep)];

end


% close all
RatioBinned = [];
RatioSEM = [];

for n = 2:20
    ii = find(nAllSubj == n);
    ratio = log10(CoupledAllSubj(ii)./unCoupledAllSubj(ii));
%     ratio = CoupledAllSubj(ii)./unCoupledAllSubj(ii);
    
    RatioBinned = [RatioBinned, mean(ratio)];
    RatioSEM = [RatioSEM, std(ratio)];%/sqrt(length(ii))];

        
    
end

h4  =  figure('Position',[288 873 632/2 476/2]);
[pl pa] = boundedline(2:20,RatioBinned,RatioSEM, '-s'); hold on;
pl.LineWidth = 1.0;
pl.Color = [0.5 0.5 0.5]; %clr(subj,:); %'k';
pl.MarkerFaceColor = [0.9 0.9 0.9];
pa.FaceColor = [0.5 0.5 0.5]; %'k'; %clr(subj,:); %'k';
pa.FaceAlpha = 0.3;
% save(fullfile(fdir,'E1_CoRipleCoSpike.mat'), 'coSpikeMatrixCoupled', 'coSpikeMatrixUncoupled', '-v7.3')
xlabel('# co Occuring ripples with spikes')
ylabel('observed prob / chance prob')

ylim([0 10])
ax = gca;
% ytick = ax.YTick;
% ytick = num2cell(ytick);
% ax.YTickLabel = ytick;

fig = gcf;
fig.Color = 'w';

savepdf(gcf, fullfile(exportDirec, sprintf('%s_coFireProbChanceProb.pdf', state)))




















