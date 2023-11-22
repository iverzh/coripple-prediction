
subjects = {'B1', 'E1', 'E2'};
state = 'wake';
matExportFolder = '';
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
for s = 1:length(subjects)
    subject = subjects{s};
    stageFolder ='';
    load(fullfile(stageFolder,'')) %load stage  mask

    if s == 1
        units = load('');
        load(fullfile(matExportFolder,'')); %load ripple detection
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
            units = load('');
            load(fullfile(matExportFolder,''));
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
        units = load('');
        load(fullfile(matExportFolder,''));
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
exportDirec = '../../figures';

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


%% unit Phase locking to ripple
close all
subject = 'B1';
state = 'NREM';
stageFolder = sprintf('');

units = load('');
load(fullfile(stageFolder,'')) %Stage Mask

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

        uT = uT(rippMask(ch,uT));
        unitPhase = RBphaseAll(ch, uT);
        phase{ui} = unitPhase;

        s = sum(unitPhase<pi/2 & unitPhase>-pi/2);
        n = s + sum(unitPhase>pi/2 | unitPhase<-pi/2);
        p = myBinomTest(s,n,0.5);

        pBinom(ui) = p;
        
    end

 

end
FRrip(FRrip <= 0) = nan;    
yy = (FRrip - FRbase) ./ FRbase;

    

ti = strcmp(type, 'pyr');
sum(pBinom(ti) < .05) / sum(ti)

ti = strcmp(type, 'int');
sum(pBinom(ti) < .05) / sum(ti)







