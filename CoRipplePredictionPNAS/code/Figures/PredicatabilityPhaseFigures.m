
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

switch state
    case 'wake'
        printFig = 3;
    case 'NREM'
        printFig = length(subjects);
end
%% Predictability vs Distance
close all
tic
errorCount = 0;
fO = [];
fOBinned = [];
fO_per_cell = [];
pO_A = [];
pO_B = [];
peak_all_subj = [];
trough_all_subj = [];
plv_A_AllSubjects = [];
plv_B_AllSubjects = [];
pp_AllSubjects = [];
for subj = 1% [1 2 3 4]% [1 2 3 4] %[1 2 3 4] %:length(subjects)
    subject = subjects{subj};
   
    try
        if contains(subject, 'E')
            arrayConfig = 'N/A';

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
  

    
    for tThresh = 150
        N = 5;
        Nplot = 5;
        bins = linspace(-pi,pi,N+1);
        ppAll = [];
        plv_A_All = [];
        plv_B_All = [];
        pAll = [];

        fO_subj = [];
        pO_A_subj = [];
        pO_B_subj = [];
    
        fO_all = [];
        pO_A_all = [];
        pO_B_all = [];
        timesR_all = [];
    
        mu_all = [];
        peak_all = [];
        trough_all = [];
        fO_cell_all = [];
        
        kl = [];
        for iiA = 1:size(rippStats.PLV_A,1)
            for iiB = 1:size(rippStats.PLV_A,2)
                if ~isempty(rippStats.FilterOut{iiA,iiB}) && ~isempty(nullStats.FilterOut{iiA,iiB})
                    
                    rippPhaseB = rippStats.PhaseB{iiA,iiB};
                    rippPhaseA = rippStats.rippPhaseUnitA{iiA,iiB};
                    s = sum(rippPhaseB<pi/2 & rippPhaseB>-pi/2);
                    pBinom = myBinomTest(s,length(rippPhaseB),0.5);
                    if isfield(rippStats,'timesWin')
                        timesR = rippStats.timesWin{iiA,iiB};
                    else
                        timesR = rippStats.times{iiA,iiB}';
                        timesR = 751 - timesR + 1;
                    end
                    filterOut = rippStats.FilterOut{iiA,iiB};
   
                    timesR(timesR > win) = [];
    
                    ppr = mean(rippStats.FilterOut{iiA,iiB}, 'omitnan');
                    ppn = mean(nullStats.FilterOut{iiA,iiB}, 'omitnan');
                    unitsAll = [rippStats.unitsA; rippStats.unitsB];
                    APar = rippStats.APar{iiA,iiB};
                    if  ppr > 0 && rippStats.Sig(iiA,iiB) > 0
                        plvA = rippStats.PLV_A(iiA,iiB);
                        phases = rippStats.PhaseB{iiA,iiB};
                        filterOut = filterOut(~isnan(phases));
                        phases = phases(~isnan(phases));
                        plvB = PLV(phases, zeros(1,length(phases)));
                        
                        ppAll = [ppAll ppr];

                        plv_A_All = [plv_A_All plvA];
                        plv_B_All = [plv_B_All plvB];
                                                    
                        pBinned = [];
                        for b = 1:length(bins)-1
                            ii = phases > bins(b) & phases < bins(b+1);
                            pBinned = [pBinned mean(filterOut(ii))];
                        
                        
                        end
            
                        pBinned = pBinned/max(pBinned);
      
                        fOBinned = [fOBinned; pBinned];

                       
                    end
                    if strcmp(subject,'B1') && strcmp(state,'wake')
                        stageKeep = find(rippStats.APar{iiA, iiB} > 3e7);
                        stageKeep = ismember(rippStats.trials{iiA, iiB}(rippStats.times{iiA, iiB} > 601), stageKeep)';

                    else
                        stageKeep = true(size(rippStats.FilterOut{iiA,iiB}));
                    end
                    interactionType = rippStats.CellType{iiA,iiB};
                    fO_subj = [fO_subj, mean(rippStats.FilterOut{iiA,iiB}(stageKeep), 'omitnan')];
                    pTemp = circ_mean(rippStats.rippPhaseUnitA{iiA,iiB}(stageKeep)');
                    pO_A_subj = [pO_A_subj, pTemp];
                    pTemp = circ_mean(rippStats.PhaseB{iiA,iiB}(stageKeep)');
                    pO_B_subj = [pO_B_subj, pTemp];

                    fO_all = [fO_all rippStats.FilterOut{iiA,iiB}(stageKeep)];
                    fO_per_cell = [fO_per_cell mean(rippStats.FilterOut{iiA,iiB}(stageKeep))];
                    pO_A_all = [pO_A_all rippStats.rippPhaseUnitA{iiA,iiB}(stageKeep)];
                    pO_B_all = [pO_B_all rippStats.PhaseB{iiA,iiB}(stageKeep)];
                    timesR_all = [timesR_all timesR(timesR <= win & (stageKeep))];
                        
                      

                    d = 3;
                    
                    fO_cell =  rippStats.FilterOut{iiA,iiB};
                    pO_A_cell = rippStats.rippPhaseUnitA{iiA,iiB};
                    pO_B_cell = rippStats.PhaseB{iiA,iiB};
                    iiB_peak = pO_B_cell > (d-1)*pi/d |  pO_B_cell < -(d-1)*pi/d;
                    iiA_peak = pO_A_cell > (d-1)*pi/d |  pO_A_cell < -(d-1)*pi/d;
                    ii_times = timesR <= tThresh; %44; % & timesR_all > (t-10);
                    iiA_trough = pO_A_cell < pi/d &  pO_A_cell > -pi/d;
                    iiB_trough = pO_B_cell < pi/d &  pO_B_cell > -pi/d;
                    
                    nPeak = sum(iiA_peak & iiB_peak & ii_times & stageKeep);
                    nTrough = sum(iiA_trough & iiB_trough & ii_times & stageKeep);
                    if nPeak > 2 && nTrough > 2
                        fO_peak = fO_cell(iiA_peak & iiB_peak & ii_times & stageKeep);
                        fO_trough = fO_cell(iiA_trough & iiB_trough & ii_times & stageKeep);
                        peak_all = [peak_all mean(fO_peak)];
                        trough_all = [trough_all mean(fO_trough)];
                        fO_cell_all = [fO_cell_all mean(fO_cell)];
                    end
                 
                        
                    
                end
            end
        end

    end
    
    pBinnedAllA = [];
    SEMall = [];
    
    for b = 1:length(bins)-1
        ii = pO_A_subj > bins(b) & pO_A_subj < bins(b+1);
        pBinnedAllA = [pBinnedAllA mean(fO_subj(ii), 'omitnan')];
        SEMall = [SEMall std(fO_subj(ii), 'omitnan')/sqrt(sum(ii))];
                        
    end

    pBinnedAllB = [];
    SEMall = [];
    
    for b = 1:length(bins)-1
        ii = pO_B_subj > bins(b) & pO_B_subj < bins(b+1);
        pBinnedAllB = [pBinnedAllB mean(fO_subj(ii), 'omitnan')];
        SEMall = [SEMall std(fO_subj(ii), 'omitnan')/sqrt(sum(ii))];
                        
    end

    figure(200); 
    subplot(2,2,subj)
    plot(pBinnedAllA, 'g-'); hold on;
    plot(pBinnedAllB, 'r-'); hold on;
    title(length(fO_subj))
    
    fO = [fO, fO_subj];
    pO_A = [pO_A, pO_A_subj];
    pO_B = [pO_B, pO_B_subj];
    
    if printFig == subj
        pO_A(isnan(fO)) = []; pO_B(isnan(fO)) = []; fO(isnan(fO)) = [];

        [sFit, modulationOrig, modulation, fvalOrig, fval] = fitsinecurve(pO_A,fO,N);
       
        figure(3); 
        fig = gcf;
        fig.Position = [425 578 562 408];
        hs = histogram(abs(modulation),100,'Normalization','probability'); hold on;
        hs.FaceColor = [0 108 0]/255; 
    
        v = vline(abs(modulationOrig)); hold on
        v.Color = 'r';
        v.LineStyle = '-';
        v.LineWidth = 2;
        
        pModulation = findPercentile(modulation', abs(modulationOrig),'pval');
        txt = sprintf('p = %.4f',pModulation);
        yl = get(gca,'ylim');
        xl = get(gca,'xlim');
        t = text(0.7*xl(2),0.9*yl(2),txt); hold on;
    
        fig.Color = 'w';
        fname = sprintf('%s_MI_UnitA_SmoothWindow_%i_PSTHwindow_%i.pdf',state,smoothingWindow, win);
        savepdf(gcf, fullfile(exportDirec,fname));
    
    
        bins = linspace(-pi,pi,Nplot+1);
    
        pBinnedAll = [];
        SEMall = [];
        
        for b = 1:length(bins)-1
            ii = pO_A > bins(b) & pO_A < bins(b+1) & isfinite(fO);
            pBinnedAll = [pBinnedAll mean(fO(ii), 'omitnan')];
            SEMall = [SEMall std(fO(ii), 'omitnan')/sqrt(sum(ii))];
                            
        end
    
        pBinnedAll = pBinnedAll/max(pBinnedAll(:));

    
        figure(1);
        
        binsHist = linspace(-pi,pi,20+1);

        binCenters = smoothdata(bins,'movmean',2);
        binCenters(1) = [];

        
        [ba, bf] = boundedline(binCenters, pBinnedAll, SEMall, '-'); hold on;
        ba.Color = [0 108 0]/255;
        bf.FaceColor = [0 108 0]/255;
        ba.LineWidth = 2;
        bf.FaceAlpha = 0.3;
        if pModulation < 0.05
    
            xp = linspace(-pi,pi);
            per = 2*pi;                     % Estimate period
            fit = @(b,x)  b(1).*(sin(2*pi*x./per + 2*pi/b(2))) + b(3);    % Function to fit
            pl = plot(xp,fit(sFit,xp), '--'); hold on;
            pl.Color = [0 108 0]/255;

            pl.LineWidth = 2;
    
            prefPhase = xp(fit(sFit,xp) == max(fit(sFit,xp)));
            
    
    
            
        end
        minLim = floor(min(pBinnedAll-SEMall)*10)/10;
        maxLim = 1.05; 
        ylim([minLim maxLim])
        
        xlim([-pi pi])  
        if strcmp(state, 'NREM')
            ylim([0.75 1])
        else
            ylim([0.2 1])
        end
        box off
        ax.YAxis(1).Color = [0 108 0]/255;
    
        xlabel('target spike ripple phase [rad]')
        ylabel('predicatbility (normalized)')

        
    

        fig = gcf;
        fig.Color = 'w';
        fig.Position =  [951 276 439 311];

        fname = sprintf('%s_PredictabilityPhaseTuning_UnitA_SmoothWindow_%i_PSTHwindow_%i_noHist.pdf',state,smoothingWindow, win);
        savepdf(gcf, fullfile(exportDirec,fname));

        [sFit, modulationOrig, modulation, fvalOrig, fval] = fitsinecurve(pO_B,fO,N);
    
        figure(5); 
     
        hs = histogram(abs(modulation),100,'Normalization','probability'); hold on;
        hs.FaceColor = [238 34 18]/255;
    
        v = vline(abs(modulationOrig)); hold on
        v.Color = 'r';
        v.LineStyle = '-';
        v.LineWidth = 2;
        
        pModulation = findPercentile(modulation', abs(modulationOrig),'pval');
        txt = sprintf('p = %.4f',pModulation);
        yl = get(gca,'ylim');
        xl = get(gca,'xlim');
        t = text(0.7*xl(2),0.9*yl(2),txt); hold on;
    
        fig.Color = 'w';
        fname = sprintf('%s_MI_UnitB_SmoothWindow_%i_PSTHwindow_%i.pdf',state,smoothingWindow, win);
        savepdf(gcf, fullfile(exportDirec,fname));
    
    
        bins = linspace(-pi,pi,Nplot+1);
    
        pBinnedAll = [];
        SEMall = [];
        
        for b = 1:length(bins)-1
            ii = pO_B > bins(b) & pO_B < bins(b+1) & isfinite(fO);
            pBinnedAll = [pBinnedAll mean(fO(ii), 'omitnan')];
            SEMall = [SEMall std(fO(ii), 'omitnan')/sqrt(sum(ii))];
                            
        end
    
        pBinnedAll = pBinnedAll/max(pBinnedAll(:));
    
        
        figure(6);
                
        yyaxis right
        binsHist = linspace(-pi,pi,20+1);
    
        hst = histogram(pO_B, binsHist, 'Normalization', 'probability'); hold on;
        colr = 1.8*[238 34 18]/255;
        colr(colr > 1) = 1;
        hst.FaceColor = [0.8 0.8 0.8];
        hst.EdgeColor = [238 34 18]/255;
        ax = gca;
        ax.YAxis(2).Color = 'k';
        ylabel('probability')
    
    
        yyaxis left
        binCenters = smoothdata(bins,'movmean',2);
        binCenters(1) = [];
        
        [ba, bf] = boundedline(binCenters, pBinnedAll, SEMall, '-'); hold on;
        ba.Color = [238 34 18]/255;
        bf.FaceColor = [238 34 18]/255;
        ba.LineWidth = 2;
        bf.FaceAlpha = 0.3;
        if pModulation < 0.05
    
            xp = linspace(-pi,pi);
            per = 2*pi;                     % Estimate period
            fit = @(b,x)  b(1).*(sin(2*pi*x./per + 2*pi/b(2))) + b(3);    % Function to fit
            pl = plot(xp,fit(sFit,xp), '--'); hold on;
            pl.Color = [238 34 18]/255;
            pl.LineWidth = 2;
    
            prefPhase = xp(fit(sFit,xp) == max(fit(sFit,xp)));
            
    
    
            
        end
        minLim = floor(min(pBinnedAll-SEMall)*10)/10;
        maxLim = 1.05; 
        ylim([minLim maxLim])
        
        xlim([-pi pi])  
        if strcmp(state, 'NREM')
            ylim([0.75 1])
        else
            ylim([0.2 1])
        end       
        box off
   
        ax = gca;
        ax.YAxis(1).Color = [238 34 18]/255;
    
        xlabel('driver spike ripple phase [rad]')
        ylabel('predicatbility (normalized)')

        fig = gcf;
        fig.Color = 'w';
        fig.Position =  [951 276 439 311];

        fname = sprintf('%s_PredictabilityPhaseTuning_UnitB_SmoothWindow_%i_PSTHwindow_%i.pdf',state,smoothingWindow, win);
        savepdf(gcf, fullfile(exportDirec,fname));
    end
   
    
    figure(2)
    fprintf('%s\n', subject)

    [RHO,PVAL] = corr([plv_A_All+plv_B_All]',ppAll','Type','Pearson');
    string = sprintf('summed PLV\n r = %.2f\np = %.3f\n', RHO, PVAL);
    fprintf(string)

    [RHO,PVAL] = corr(plv_A_All',ppAll','Type','Pearson');
    string = sprintf('target PLV\n r = %.2f\np = %.3f\n', RHO, PVAL);
    fprintf(string)

    [RHO,PVAL] = corr(plv_B_All',ppAll','Type','Pearson');
    string = sprintf('driver PLV\n r = %.2f\np = %.3f\n', RHO, PVAL);
    fprintf(string)
    [p, S, CI ]= polyfitconfidence(plv_A_All,ppAll);
    x1 = linspace(0,max(plv_A_All));
    [y1,delta] = polyval(p,x1,S);

    plv_A_AllSubjects = [plv_A_AllSubjects, plv_A_All];
    plv_B_AllSubjects = [plv_B_AllSubjects, plv_B_All];
    pp_AllSubjects = [pp_AllSubjects, ppAll];
    peak_all_subj = [peak_all_subj, peak_all];
    trough_all_subj = [trough_all_subj, trough_all];
    
    pl = plot(x1, y1, '-'); hold on;
    pl.LineWidth = 3;
    pl.Color = clr(subj,:);

    pl = plot(x1, CI(1,:), '--'); hold on;
    pl.LineWidth = 1.5;
    pl.Color = clr(subj,:);

    pl = plot(x1, CI(2,:), '--'); hold on;
    pl.LineWidth = 1.5;
    pl.Color = clr(subj,:);

    
    xlabel('ripple target-spike PLV')
    ylabel('predicatbility')
    xlim([0.0 0.8])
    ylim([0.0 0.8])

    
    fig = gcf;
    fig.Position = [26 -141 398 490];
    if subj == printFig
        
        box off
        ax = gca;
        fig.Color = 'w';

        fname = sprintf('%s_PredictabilityPLVcorr_SmoothWindow_%i_PSTHwindow_%i.pdf',state,smoothingWindow, win);
        savepdf(gcf, fullfile(exportDirec,fname));
    end

    if any([strcmp(subject,'E1-') strcmp(arrayConfig,'lateral-medial-')])

        figure(50); 
        fig = gcf;
        fig.Position = [867 687 310 390];
        dscatter(plv_A_All',ppAll'); hold on;
        colormap(viridis)

        pl = plot(x1, y1, '-'); hold on;
        pl.LineWidth = 3;
        pl.Color = 'k';
        xlabel('ripple driver-spike PLV')
        ylabel('predicatbility')

        fig.Color = 'w';

        ylim([0 1.5])

        fname = sprintf('%s_%s_PredictabilityPLVcorr_SmoothWindow_%i_PSTHwindow_%i.pdf',subject, state,smoothingWindow, win);
        savepdf(gcf, fullfile(exportDirec,fname));

        clf 

    end
    

    figure(20);
    bplt = bar(subj, mean((peak_all-trough_all),'omitnan')/mean(fO_cell_all,'omitnan') ); hold on
    bplt.FaceColor = clr(subj,:);

    effect = computeCohen_d(peak_all,trough_all);

    [h, p, ci, stats] = ttest(peak_all,trough_all, 'tail', 'right');

    titl = sprintf('d: %.2f\n p: %.3f', effect, p);
    if subj == printFig
        box off 

        ax= gca;
        ax.XTick = 1:length(subjects);
        
        fig = gcf;
        fig.Position =  [435 257 295 464*0.5];

        fname = sprintf('%s_TroughPeakPredict_SmoothWindow_%i_PSTHwindow_%i.pdf',state,smoothingWindow, win);
        savepdf(gcf, fullfile(exportDirec,fname));
    end
    


end

fig = gcf;
fig.Color = 'w';
ylim([0 0.5])

[RHO,PVAL] = corr([plv_A_AllSubjects+plv_B_AllSubjects]',pp_AllSubjects','Type','Pearson');
string = sprintf('ALL SUBJECTS   \n summed PLV\n r = %.2f\np = %.3f\n', RHO, PVAL);
fprintf(string)

[RHO,PVAL] = corr(plv_A_AllSubjects',pp_AllSubjects','Type','Pearson');
string = sprintf('target PLV\n r = %.2f\np = %.3f\n', RHO, PVAL);
fprintf(string)

[RHO,PVAL] = corr(plv_B_AllSubjects',pp_AllSubjects','Type','Pearson');
string = sprintf('driver PLV\n r = %.2f\np = %.3f\n', RHO, PVAL);
fprintf(string)


[h, p, ci, stats] = ttest(peak_all_subj, trough_all_subj, 'tail', 'right');

effect = computeCohen_d(peak_all_subj,trough_all_subj);



%%

% Then we bootstrap the data
b = plv_B_AllSubjects';
c = [plv_A_AllSubjects+plv_B_AllSubjects]';
a = pp_AllSubjects';
d = plv_A_AllSubjects';

Nb = 1000;
Np = length(a);
bootcorr1 = zeros(Nb,1);
bootcorr2 = zeros(Nb,1);
bootcorr3 = zeros(Nb,1);
for B = 1:Nb

    bootsample = randi(Np,1,Np);
    bootcorr1(B) = corr(a(bootsample),b(bootsample),'Type','Pearson');
    bootcorr2(B) = corr(a(bootsample),d(bootsample),'Type','Pearson');
    bootcorr3(B) = corr(a(bootsample),c(bootsample),'Type','Pearson');

end

figure; 
h = histogram(bootcorr3 - bootcorr2, -0.05:0.005:0.2); hold on;
h.FaceColor = [0 108 0]/255;          
h = histogram(bootcorr3 - bootcorr1, -0.05:0.005:0.2); hold on;
h.FaceColor = [238 34 18]/255;
vline(0)
xlim([-0.05 0.2]);
box off;
fig = gcf;
fig.Color = 'w';

CI = quantile(bootcorr3 - bootcorr1, [0.025 0.5 0.975]);
fprintf('bootstap for driver vs combined %.2f [%.2f %.2f]\n', CI(2), CI(1), CI(3))
CI = quantile(bootcorr3 - bootcorr2, [0.025 0.5 0.975]);
fprintf('bootstap for target vs combined %.2f [%.2f %.2f]\n', CI(2), CI(1), CI(3))


fname = sprintf('%s_BootstappedPLVcorr_SmoothWindow_%i_PSTHwindow_%i.pdf',state,smoothingWindow, win);
savepdf(gcf, fullfile(exportDirec,fname));






%%
binsB = quantile(plv_B_AllSubjects, 9);
binsB(2:end+1) = binsB; binsB(1) = min(plv_B_AllSubjects); binsB(end+1) = 1;
binsA = quantile(plv_B_AllSubjects, 9);
binsA(2:end+1) = binsA; binsA(1) = min(plv_A_AllSubjects); binsA(end+1) = 1;

predictMat = nan(length(binsA)-1);
for Ab = 1:length(binsA)-1
    for Bb = 1:length(binsB)-1
        iiA = plv_A_AllSubjects >= binsA(Ab) & plv_A_AllSubjects < binsA(Ab+1);
        iiB = plv_B_AllSubjects >= binsB(Bb) & plv_B_AllSubjects < binsB(Bb + 1);
        ii = iiA & iiB; 
        predictMat(Ab, Bb) = mean(pp_AllSubjects(ii), 'omitnan');   
    end
    
end
figure; imagesc(flipud(predictMat)); colormap(inferno); colorbar
fname = sprintf('%s_PLVpredictMatrix_SmoothWindow_%i_PSTHwindow_%i.pdf',state,smoothingWindow, win);
savepdf(gcf, fullfile(exportDirec,fname));










