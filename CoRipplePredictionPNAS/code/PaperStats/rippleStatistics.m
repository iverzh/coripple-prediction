










subjects = {'T11', 'MG29', 'MG63'};
state = 'wake';
matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
%
density = [];
amp = [];
duration = [];
freq = [];
PYripSpikes = [];
INripSpikes = [];
PYfrRip = [];
INfrRip = [];
PYfrNull = [];
INfrNull = [];
for s = 1 %2:length(subjects)
    subject = subjects{s};
    stageFolder = sprintf('/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/%s/data_1kHz', subject);
    load(fullfile(stageFolder, sprintf('%s_stageMask_%s.mat', subject, state)))

    if s == 1
        units = LoadSpikeTimes(subject,'CellExplorer','medial');
        load(fullfile(matExportFolder,[subject,'_ripple_stats_wake_NC_wakeonly_medial.mat']));
%         load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMonly_medial.mat']));
        channels = unique(cell2mat(units(:,1)));
        for ii = 1:length(channels)
            ch = channels(ii);
            amp = [amp mean(rippleStats.rippleAmp{ch})];
            duration = [duration mean(rippleStats.duration{ch})];
            freq = [freq mean(rippleStats.oscFreq{ch})];
            
            locs = rippleStats.locs{ch};
            locs = locs(StageMask_1kHz(round(locs)));
            density = [density length(locs)/sum(StageMask_1kHz) * 1e3 * 60];

            unitsOnCh = cellfun(@(X) X == ch, units(:,1));
            types = cellfun(@(X) strcmp(X,'pyr'), units(:,3));
            AP_ch = round(cell2mat(units(unitsOnCh & types,2)) * 1e3);
            nSpk = [];
            for r = 1:length(rippleStats.locs{ch})
                window = rippleStats.window{ch}(r,:);
                nSpk = [nSpk sum(AP_ch >= window(1) & AP_ch <= window(2))];
            end

            PYripSpikes = [PYripSpikes mean(nSpk(nSpk > 0))];
            FR = sum(nSpk) / sum(rippleStats.duration{ch}) * 1e3;
            PYfrRip =[PYfrRip FR];
            AP_ch = AP_ch(StageMask_1kHz(AP_ch));
            FR = (length(AP_ch) - sum(nSpk)) / (sum(StageMask_1kHz) - sum(rippleStats.duration{ch})) * 1e3;
            PYfrNull =[PYfrNull FR];

            types = cellfun(@(X) strcmp(X,'int'), units(:,3));
            AP_ch = round(cell2mat(units(unitsOnCh & types,2)) * 1e3);
            nSpk = [];
            for r = 1:length(rippleStats.locs{ch})
                window = rippleStats.window{ch}(r,:);
                nSpk = [nSpk sum(AP_ch >= window(1) & AP_ch <= window(2))];
            end

            INripSpikes = [INripSpikes mean(nSpk(nSpk > 0))];
            FR = sum(nSpk) / sum(rippleStats.duration{ch}) * 1e3;
            INfrRip =[INfrRip FR];
            AP_ch = AP_ch(StageMask_1kHz(AP_ch));
            FR = (length(AP_ch) - sum(nSpk)) / (sum(StageMask_1kHz) - sum(rippleStats.duration{ch})) * 1e3;
            INfrNull =[INfrNull FR];
        end
        if strcmp(state, 'NREM')
            units = LoadSpikeTimes(subject,'CellExplorer','lateral');
            load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMonly_lateral.mat']));
            channels = unique(cell2mat(units(:,1)));
            for ii = 1:length(channels)
                ch = channels(ii);
                amp = [amp mean(rippleStats.rippleAmp{ch})];
                duration = [duration mean(rippleStats.duration{ch})];
                freq = [freq mean(rippleStats.oscFreq{ch})];
    
                locs = rippleStats.locs{ch};
                locs = locs(StageMask_1kHz(round(locs)));
                density = [density length(locs)/sum(StageMask_1kHz) * 1e3 * 60];
    
                unitsOnCh = cellfun(@(X) X == ch, units(:,1));
                types = cellfun(@(X) strcmp(X,'pyr'), units(:,3));
                AP_ch = round(cell2mat(units(unitsOnCh & types,2)) * 1e3);
                nSpk = [];
                for r = 1:length(rippleStats.locs{ch})
                    window = rippleStats.window{ch}(r,:);
                    nSpk = [nSpk sum(AP_ch >= window(1) & AP_ch <= window(2))];
                end
    
                PYripSpikes = [PYripSpikes mean(nSpk(nSpk > 0))];
                FR = sum(nSpk) / sum(rippleStats.duration{ch}) * 1e3;
                PYfrRip =[PYfrRip FR];
                AP_ch = AP_ch(StageMask_1kHz(AP_ch));
                FR = (length(AP_ch) - sum(nSpk)) / (sum(StageMask_1kHz) - sum(rippleStats.duration{ch})) * 1e3;
                PYfrNull =[PYfrNull FR];
    
                types = cellfun(@(X) strcmp(X,'int'), units(:,3));
                AP_ch = round(cell2mat(units(unitsOnCh & types,2)) * 1e3);
                nSpk = [];
                for r = 1:length(rippleStats.locs{ch})
                    window = rippleStats.window{ch}(r,:);
                    nSpk = [nSpk sum(AP_ch >= window(1) & AP_ch <= window(2))];
                end
    
                INripSpikes = [INripSpikes mean(nSpk(nSpk > 0))];
                FR = sum(nSpk) / sum(rippleStats.duration{ch}) * 1e3;
                INfrRip =[INfrRip FR];
                AP_ch = AP_ch(StageMask_1kHz(AP_ch));
                FR = (length(AP_ch) - sum(nSpk)) / (sum(StageMask_1kHz) - sum(rippleStats.duration{ch})) * 1e3;
                INfrNull =[INfrNull FR];
     
            end
        end
    else
        units = LoadSpikeTimes(subject,'CellExplorer');
        load(fullfile(matExportFolder,[subject,'_ripple_stats_wake_NC_wakeOnly.mat']));
%         load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMOnly.mat']));
        channels = unique(cell2mat(units(:,1)));
        for ii = 1:length(channels)
            ch = channels(ii);
            amp = [amp mean(rippleStats.rippleAmp{ch})];
            duration = [duration mean(rippleStats.duration{ch})];
            freq = [freq mean(rippleStats.oscFreq{ch})];

            locs = rippleStats.locs{ch};
            locs = locs(StageMask_1kHz(round(locs)));
            density = [density length(locs)/sum(StageMask_1kHz) * 1e3 * 60];

            unitsOnCh = cellfun(@(X) X == ch, units(:,1));
            types = cellfun(@(X) strcmp(X,'pyr'), units(:,3));
            AP_ch = round(cell2mat(units(unitsOnCh & types,2)) * 1e3);
            nSpk = [];
            for r = 1:length(rippleStats.locs{ch})
                window = rippleStats.window{ch}(r,:);
                nSpk = [nSpk sum(AP_ch >= window(1) & AP_ch <= window(2))];
            end

            PYripSpikes = [PYripSpikes mean(nSpk(nSpk > 0))];
            FR = sum(nSpk) / sum(rippleStats.duration{ch}) * 1e3;
            PYfrRip =[PYfrRip FR];
            AP_ch = AP_ch(StageMask_1kHz(AP_ch));
            FR = (length(AP_ch) - sum(nSpk)) / (sum(StageMask_1kHz) - sum(rippleStats.duration{ch})) * 1e3;
            PYfrNull =[PYfrNull FR];
            
            types = cellfun(@(X) strcmp(X,'int'), units(:,3));
            AP_ch = round(cell2mat(units(unitsOnCh & types,2)) * 1e3);
            nSpk = [];
            for r = 1:length(rippleStats.locs{ch})
                window = rippleStats.window{ch}(r,:);
                nSpk = [nSpk sum(AP_ch >= window(1) & AP_ch <= window(2))];
            end

            INripSpikes = [INripSpikes mean(nSpk(nSpk > 0))];
            INripSpikes = [INripSpikes mean(nSpk(nSpk > 0))];
            FR = sum(nSpk) / sum(rippleStats.duration{ch}) * 1e3;
            INfrRip =[INfrRip FR];
            AP_ch = AP_ch(StageMask_1kHz(AP_ch));
            FR = (length(AP_ch) - sum(nSpk)) / (sum(StageMask_1kHz) - sum(rippleStats.duration{ch})) * 1e3;
            INfrNull =[INfrNull FR];
 
        end
    end
end

fprintf('PY ripp: %0.2fHz %s %0.2fHz   control: %0.2fHz %s %0.2fHz\n', mean(PYfrRip(PYfrRip~=0)), char(177), std(PYfrRip(PYfrRip~=0)),  mean(PYfrNull(PYfrNull~=0)), char(177), std(PYfrNull(PYfrNull~=0)))
fprintf('IN ripp: %0.2fHz %s %0.2fHz   control: %0.2fHz %s %0.2fHz\n', mean(INfrRip(INfrRip~=0)), char(177), std(INfrRip(INfrRip~=0)),  mean(INfrNull(INfrNull~=0)), char(177), std(INfrNull(INfrNull~=0)))


[P,H] = signrank(PYfrRip(PYfrRip~=0), PYfrNull(PYfrRip~=0), 'tail','both')
[P,H] = signrank(INfrRip(INfrRip~=0), INfrNull(INfrRip~=0), 'tail','both')

densityM = density;
ampM = amp;
durationM = duration;
freqM = freq;
%% Compare ripple charateristics from M1 and temporal lobe
exportDirec = '/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/SpikePredict/figures/CoRippleFigures';

figure('Position', [435 368 241 413]);
M = [densityM,densityT];
cgroup = [ones(1,length(densityM)), 2*ones(1,length(densityT))];
bp = boxchart(M', 'notch', 'off', 'GroupByColor', cgroup);

bp(1).MarkerStyle = 'none';
bp(2).MarkerStyle = 'none';

bp(1).BoxFaceColor = mean(clr(1:2,:));
bp(2).BoxFaceColor = mean(clr(3:4,:));

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirec, sprintf('density_%s.pdf', state)))

figure('Position', [435 368 241 413]);
M = [ampM,ampT];
cgroup = [ones(1,length(densityM)), 2*ones(1,length(densityT))];
bp = boxchart(M', 'notch', 'off', 'GroupByColor', cgroup);

bp(1).MarkerStyle = 'none';
bp(2).MarkerStyle = 'none';

bp(1).BoxFaceColor = mean(clr(1:2,:));
bp(2).BoxFaceColor = mean(clr(3:4,:));

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirec, sprintf('amplitude_%s.pdf', state)))

figure('Position', [435 368 241 413]);
M = [durationM,durationT];
cgroup = [ones(1,length(densityM)), 2*ones(1,length(densityT))];
bp = boxchart(M', 'notch', 'off', 'GroupByColor', cgroup);

bp(1).MarkerStyle = 'none';
bp(2).MarkerStyle = 'none';

bp(1).BoxFaceColor = mean(clr(1:2,:));
bp(2).BoxFaceColor = mean(clr(3:4,:));

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirec, sprintf('duration_%s.pdf', state)))

figure('Position', [435 368 241 413]);
M = [freqM,freqT];
cgroup = [ones(1,length(densityM)), 2*ones(1,length(densityT))];
bp = boxchart(M', 'notch', 'off', 'GroupByColor', cgroup);

bp(1).MarkerStyle = 'none';
bp(2).MarkerStyle = 'none';

bp(1).BoxFaceColor = mean(clr(1:2,:));
bp(2).BoxFaceColor = mean(clr(3:4,:));

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirec, sprintf('frequency_%s.pdf', state)))


%% unit Phase locking 
close all
subject = 'T11';
state = 'NREM';
stageFolder = sprintf('/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/%s/data_1kHz', subject);

units = LoadSpikeTimes(subject,'CellExplorer','lateral');
load(fullfile(stageFolder, sprintf('%s_stageMask_%s.mat', subject, state)))
% load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMOnly.mat']));
% load(fullfile(matExportFolder,[subject,'_ripple_stats_wake_NC_wakeOnly.mat']));
% load(fullfile(matExportFolder,[subject,'_ripple_stats_wake_NC_wakeonly_medial.mat']));
% load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMonly_medial.mat']));
load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMonly_lateral.mat']));

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

rippMask = logical(rippMask);
rippMask(:, ~StageMask_1kHz) = false;
type = [];
phase = [];
pBinom = [];
% FRrip = [];
% FRbase = [];
for ui = 1:length(units)
    if any(strcmp(units{ui, 3},  {'pyr', 'int'}))
        uT = round(units{ui,2} * rippleStats.fs);
        ch = units{ui,1};
        type{ui} = units{ui,3};
        FRrip = [FRrip (sum(rippMask(ch,uT)) / sum(rippMask(ch,:))) * 1e3];
        FRbase = [FRbase (sum(StageMask_1kHz(uT)) / sum(StageMask_1kHz)) * 1e3];

%         uT = uT(rippMask(ch,uT));
%         unitPhase = RBphaseAll(ch, uT);
%         phase{ui} = unitPhase;
% 
%         s = sum(unitPhase<pi/2 & unitPhase>-pi/2);
%         n = s + sum(unitPhase>pi/2 | unitPhase<-pi/2);
%         p = myBinomTest(s,n,0.5);
% 
%         pBinom(ui) = p;
        
    end

 

end
FRrip(FRrip <= 0) = nan;    
yy = (FRrip - FRbase) ./ FRbase;

    

ti = strcmp(type, 'pyr');
sum(pBinom(ti) < .05) / sum(ti)

ti = strcmp(type, 'int');
sum(pBinom(ti) < .05) / sum(ti)

% save(sprintf('%s_SpikeRipplePhases_%s.mat', subject, state), 'phase', 'type', "pBinom", '-v7.3')
%%
NREM = load('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/SpikePredict/figures/PredictPhaseFigures/PLV/FRincreaseNREM.mat');
Wake = load('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/SpikePredict/figures/PredictPhaseFigures/PLV/FRincreaseWaking.mat');
close all
figure;
pl = plot(log(NREM.FRbase),NREM.yy , '.'); hold on;
% Step 2: Perform Linear Regression
X = log(NREM.FRbase(~isnan(NREM.yy)));
Y2 = NREM.yy(~isnan(NREM.yy));
coefficients = polyfit(X, Y2, 1); % Fit a straight line (degree 1)

% Step 3: Plot Data and Regression Line
x_fit = linspace(min(X),max(X),1000); % Generate x values for the regression line
y_fit = polyval(coefficients, x_fit); % Compute y values for the regression line
ln = plot(x_fit, y_fit, 'r', 'LineWidth', 3); hold on;% Plot the regression line in red
ln.Color = 0.8*pl.Color;

pl = plot(log(Wake.FRbase),Wake.yy , '.'); hold on;
X = log(Wake.FRbase(~isnan(Wake.yy)));
Y2 = Wake.yy(~isnan(Wake.yy));
coefficients = polyfit(X, Y2, 1); % Fit a straight line (degree 1)

% Step 3: Plot Data and Regression Line
x_fit = linspace(min(X),max(X),1000); % Generate x values for the regression line
y_fit = polyval(coefficients, x_fit); % Compute y values for the regression line
ln = plot(x_fit, y_fit, 'r', 'LineWidth', 3); hold on;% Plot the regression line in red
ln.Color = 0.8*pl.Color;

fig = gcf;
fig.Color = 'w';
fig.Position = [866 764 345 413];

ylabel('change in FR during ripples')
xlabel('log(FR)')
box off
[RHO,PVAL] = corr(NREM.FRbase(~isnan(NREM.yy))',NREM.yy(~isnan(NREM.yy))','Type','Pearson')
[RHO,PVAL] = corr(Wake.FRbase(~isnan(Wake.yy))',Wake.yy(~isnan(Wake.yy))','Type','Pearson')

savepdf(gcf,'/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/SpikePredict/figures/rippFRvsBaseline.pdf')

NREM = load('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/SpikePredict/figures/PredictPhaseFigures/PLV/phasedataNREM.mat');
Wake = load('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/SpikePredict/figures/PredictPhaseFigures/PLV/phasedataWaking.mat');
ii = NREM.pp_AllSubjects < 3;
X1 = log(NREM.FR_AllSubjects(ii));
Y1 = NREM.plv_B_AllSubjects(ii);
figure;
pl = plot(X1,Y1 , '.'); hold on;
% Step 2: Perform Linear Regression

coefficients = polyfit(X1, Y1, 1); % Fit a straight line (degree 1)

% Step 3: Plot Data and Regression Line
x_fit = linspace(min(X1),max(X1),1000); % Generate x values for the regression line
y_fit = polyval(coefficients, x_fit); % Compute y values for the regression line
ln = plot(x_fit, y_fit, 'r', 'LineWidth', 3); hold on;% Plot the regression line in red
ln.Color = 0.8*pl.Color;

ii = Wake.pp_AllSubjects < 3;
X2 = log(Wake.FR_AllSubjects(ii));
Y2 = Wake.plv_B_AllSubjects(ii); 
pl = plot(X2, Y2 , '.'); hold on;

coefficients = polyfit(X2, Y2, 1); % Fit a straight line (degree 1)

% Step 3: Plot Data and Regression Line
x_fit = linspace(min(X2),max(X2),1000); % Generate x values for the regression line
y_fit = polyval(coefficients, x_fit); % Compute y values for the regression line
ln = plot(x_fit, y_fit, 'r', 'LineWidth', 3); hold on;% Plot the regression line in red
ln.Color = 0.8*pl.Color;

fig = gcf;
fig.Color = 'w';
fig.Position = [866 764 345 413];
ylabel('PLVb')
xlabel('log(FR)')
box off
[RHO,PVAL] = corr(X1', Y1','Type','Pearson')
[RHO,PVAL] = corr(X2',Y2','Type','Pearson')
savepdf(gcf,'/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/SpikePredict/figures/PLVBvsBaseline.pdf')
ii = NREM.pp_AllSubjects < 3;
X1 = log(NREM.FR_AllSubjects(ii));
Y1 = NREM.pp_AllSubjects(ii);
figure;
pl = plot(X1,Y1 , '.'); hold on;
% Step 2: Perform Linear Regression

coefficients = polyfit(X1, Y1, 1); % Fit a straight line (degree 1)

% Step 3: Plot Data and Regression Line
x_fit = linspace(min(X1),max(X1),1000); % Generate x values for the regression line
y_fit = polyval(coefficients, x_fit); % Compute y values for the regression line
ln = plot(x_fit, y_fit, 'r', 'LineWidth', 3); hold on;% Plot the regression line in red
ln.Color = 0.8*pl.Color;

ii = Wake.pp_AllSubjects < 3;
X2 = log(Wake.FR_AllSubjects(ii));
Y2 = Wake.pp_AllSubjects(ii); 
pl = plot(X2, Y2 , '.'); hold on;

coefficients = polyfit(X2, Y2, 1); % Fit a straight line (degree 1)

% Step 3: Plot Data and Regression Line
x_fit = linspace(min(X2),max(X2),1000); % Generate x values for the regression line
y_fit = polyval(coefficients, x_fit); % Compute y values for the regression line
ln = plot(x_fit, y_fit, 'r', 'LineWidth', 3); hold on;% Plot the regression line in red
ln.Color = 0.8*pl.Color;

fig = gcf;
fig.Color = 'w';
fig.Position = [866 764 345 413];
ylabel('PLVb')
xlabel('log(FR)')
box off
[RHO,PVAL] = corr(X1', Y1','Type','Pearson')
[RHO,PVAL] = corr(X2',Y2','Type','Pearson')
savepdf(gcf,'/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/SpikePredict/figures/PredictvsBaseline.pdf')

%%
pBinomAll = [];
typeAll = [];
subjects = {'T11', 'MG29', 'MG63'};
state = 'wake'; 

for s = 1:2
    subject = subjects{s};
    folder = sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/SpikePredict/PredictAll/%s', subject);
    file = sprintf('%s_SpikeRipplePhases_%s.mat', subject, state);
    load(fullfile(folder,file))
    pBinomAll = [pBinomAll pBinom];
    typeAll = [typeAll type]; 

end

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pBinomAll,.05,'pdep','yes');

ti = strcmp(typeAll, 'pyr');
sum(adj_p(ti) < .05) / sum(ti)

ti = strcmp(typeAll, 'int');
sum(adj_p(ti) < .05) / sum(ti)








