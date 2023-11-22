
close all 
clear 
clc

addpath(genpath('../../code'))

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
                nCoSpikeCo = size(sp,2); 
                
                sp = noRip(:,noRip(1,:) > 0 & noRip(2,:) > 0);
                nCoSpikeNo = size(sp,2);
                
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
        b1 = errorbar(xVal, nCoSpikeCoDist, nCoSpikeCoDistSEM, 'o-'); hold on;
        b1.MarkerFaceColor = clr(subj,:);
        b1.Color = clr(subj,:);
        b1.LineWidth = 1.0;
        
        
        
        
        p = polyfit(xVal, nCoSpikeNoDist,1); %linear fit
        x1 = linspace(min(xVal),max(xVal));
        y1 = polyval(p,x1);
        
        b2 = errorbar(xVal, nCoSpikeNoDist, nCoSpikeNoDistSEM, 'o--'); hold on;
        b2.MarkerFaceColor = 0.6*clr(subj,:);
        b2.Color = 0.6*clr(subj,:);
        b2.LineWidth = 1.0;
        b2.MarkerSize = 5;
        
        
        ylabel('CoFire Event Rate [z-score]')

        fig = gcf;
        fig.Color = 'w';
        box off
        
        xlim([0 4000])
    
        ax = gca;
        ax.XTick = 0:500:4000;

    else
       
        figure('Position', [852.0000 681 149.0000 413]);
        
        zAll = zscore([nCoSpikeCoAll nCoSpikeNoAll]);
        nCoSpikeCoAll = zAll(1:length(nCoSpikeCoAll));
        nCoSpikeNoAll = zAll(length(nCoSpikeCoAll)+1:end);
        
        clr2 = brewermap(10,'Paired');
        
        
        b1 = plot(3, mean(nCoSpikeCoAll), 'o'); hold on;
        b1.MarkerFaceColor = clr(subj,:);
        b1.Color = clr(subj,:);
        b1.LineWidth = 1.0;
        
        
        
        
        
        b2 = plot(3, mean(nCoSpikeNoAll), 'o'); hold on;
        b2.MarkerFaceColor = 0.6*clr(subj,:);
        b2.Color = 0.6*clr(subj,:);
        b2.LineWidth = 1.0;
        
        
        ylabel('CoFire Event Rate [Hz]')
        
        fig = gcf;
        fig.Color = 'w';
        box off
        
        xlim([2.5 3.5])
    
        

    end
    
    


end
ylim([-0.4000 1]) %ax.YTick(end)])

ax = gca;


savepdf(gcf, fullfile(exportDirec, [subjects{subj},'CoRippleCoFire_NREM.pdf']))
keep = nNo > 0 & nCo > 0; 
FRco = nCo(keep) ./ durationAll(keep);
fprintf('co ripple FR: %0.3f %s %0.3f ', mean(FRco), char(177), computeSEM(FRco,1))
FRno = nNo(keep) ./ durationAll(keep);
fprintf('no ripple FR: %0.3f %s %0.3f\n', mean(FRno), char(177), computeSEM(FRno,1))



%% co fire avalanche plot 
subjects = {'B1', 'B1-dual', 'E1', 'E2'};
state = 'NREM';
clr = [237/255 110/255 87/255; %B1
       255/255 66/255 161/255; %B1-dual
       108/255 227/255 208/255; %E1
       129/255 213/255 83/255]; %E2
CoupledAllSubj = [];
unCoupledAllSubj = [];
nAllSubj = [];
for subj = [1 3]
    subject = subjects{subj};
    load(sprintf('../../out/cofire/%s_%s_CoFireAvalanche.mat', subject, state))
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
xlabel('# co Occuring ripples with spikes')
ylabel('observed prob / chance prob')

ylim([0 10])
ax = gca;


fig = gcf;
fig.Color = 'w';

savepdf(gcf, fullfile(exportDirec, sprintf('%s_coFireProbChanceProb.pdf', state)))




















