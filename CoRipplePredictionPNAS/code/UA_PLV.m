
addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))
addpath(genpath('/space/seh9/2/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh9/2/halgdev/projects/iverzh/ripples/code/util'))


home_dir = '/space/seh8/2/mdeh1-3/halgdev/projects/cgonzalez/Units';
chmap_file = sprintf('%s/mg49_and_mg29_neuroport_channel_map.csv',home_dir);
chMap = dlmread(chmap_file);
%%



subject = 'MG29';
minOverlap = 25; %ms
plvSmoothWin = 10; %ms
win =  500; %ms 
randWin = [-10000 -2000];
nShuff = 100;
f = sprintf('/space/seh9/2/halgdev/projects/iverzh/data/UtahArrayData/%s/data_1kHz/%s_1kHz_unitsRemoved_wake.mat', subject, subject);
load(f, 'data', 'chan_labels');
LFP = data; clear data;
matExportFolder = '/space/seh9/2/halgdev/projects/iverzh/ripples/matFiles';
load(fullfile(matExportFolder,[subject,'_ripple_stats_wake_NC_wakeOnly.mat']));
LFP(isnan(LFP)) = 0;

RB = zeros(size(LFP));
RBphase = zeros(size(LFP));
[bButter,aButter] = butter(3, [70 100]/(rippleStats.fs/2)); %set up butterworth filter 
for ch = 1:size(LFP,1)
    RB(ch,:) = filtfilt(bButter,aButter,LFP(ch,:));
    RBphase(ch,:) = angle(hilbert(RB(ch,:)));
    fprintf('filterning and calculating phase for channel %i\n', ch)

end

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

%%
nCh = size(LFP,1);
coOccur = nan(nCh,nCh);
for chA = 1:nCh
    for chB = 1:nCh
        rippMaskAB = rippMask(chA,:) & rippMask(chB,:);
        b = mask2bounds(rippMaskAB);

        coOccurRate = size(b,1) / length(rippleStats.locs{chA});

        coOccur(chA,chB) = coOccurRate;

        
    end

    fprintf('%i\n', chA)

end
%%
PLVstruct = [];
nCh = size(RBphaseAll,1);
maxPLV = nan(nCh,nCh);
for chA = 1:nCh
    for chB = 1:nCh
        clc
        fprintf('calculating plv for %i , %i \n', chA, chB)
        
        rippOverlap = rippMask(chA,:) & rippMask(chB,:);
        b = mask2bounds(rippOverlap);
        overlapDuration = b(:,2) - b(:,1);
        b = b(overlapDuration > minOverlap,:);

        coRippleCenters = round(mean(b,2));

        rippA = zeros(length(coRippleCenters), 2*win+1);
        rippB = zeros(length(coRippleCenters), 2*win+1);
        phaseLags = zeros(length(coRippleCenters), 1);

        shufA = zeros(length(coRippleCenters), 2*win+1, nShuff);
        shufB = zeros(length(coRippleCenters), 2*win+1, nShuff);
        for r = 1:length(coRippleCenters)
            cntr = coRippleCenters(r);
            if cntr < abs(randWin(1))
                continue
            end

            phaseA = RBphaseAll(chA, cntr-win : cntr+win);
            phaseB = RBphaseAll(chB, cntr-win : cntr+win);

            rippA(r,:) = phaseA;
            rippB(r,:) = phaseB;

            phaseAzoom = RBphaseAll(chA, cntr-(win/10) : cntr+(win/10));
            phaseBzoom = RBphaseAll(chB, cntr-(win/10) : cntr+(win/10));
            
            
            phaseLags(r) = circ_mean(angdiff(phaseAzoom, phaseBzoom)');


            for n = 1:nShuff
                randCntr = 0;
                while randCntr - win <= 0 || randCntr + win > rippleStats.recordingLength
                    randCntr = cntr + randi(randWin, 1, 1); % within param.rand_win of each ripple
                end
    
                phaseA = RBphaseAll(chA, randCntr-win : randCntr+win);
                phaseB = RBphaseAll(chB, randCntr-win : randCntr+win);
    
                shufA(r,:,n) = phaseA;
                shufB(r,:,n) = phaseB;
            end
        end

        plv = PLV(rippA', rippB');
        plvSmooth = smoothdata(plv, 'gaussian', plvSmoothWin);
        PLVstruct.plv{chA,chB} = plvSmooth;
        PLVstruct.phaseLags{chA,chB} = phaseLags;

        for n = 1:nShuff
            plvShuff = PLV(shufA(r,:,n)', shufB(r,:,n)');
            plvSmoothShuff = smoothdata(plvShuff, 'gaussian', plvSmoothWin);

            PLVstruct.plvShuff{chA,chB,n} = plvSmoothShuff;
        end

        maxPLV(chA,chB) = max(plvSmooth);

        
        
    end

           
end

save(sprintf('%s_PLV_CoOccur.mat', subject), "PLVstruct", "maxPLV", "coOccur", '-v7.3')

%%
close all
% plvShuff = PLV(shufA(:,:,n)', shufB(:,:,n)');
% plv = PLV(rippA', rippB');
% plvSmooth = smoothdata(plv, 'gaussian', plvSmoothWin);

% close all
% 
% figure; 
% plot(plv); hold on;
% plot(plvShuff);

badChannels = [7,16,17,33,36,43,48,50,52,54,58,72,79,87]; 
bad = false(1,96);
bad(badChannels) = true;
distAll = [];
for ch = 1:96
    d = findContactDistance(chMap, ch, 1);
    distAll = [distAll, d];
end

e = [];
figure;
for it = 1:1000
    for nF = 1:30
    %     nF = 10;
        nnmfPLV = coRipPLV;
        nnmfPLV(:, bad) = [];
        nnmfPLV(bad,:) = [];
        diagii = logical(diag(ones(1,size(nnmfPLV,2))));
        nnmfPLVo = nnmfPLV;
        nnmfPLV(diagii) = [];
        nnmfPLV = reshape(nnmfPLV, sum(~bad)-1, sum(~bad))';
        
        [W, H] = nnmf(nnmfPLV, nF);
        err = nnmfPLV - W*H;
        e(nF,it) = sum(err(:).^2);
    
    end
    
    plot(e(:,it), 'k-'); hold on;
end
pl = plot(mean(e,2), 'r-'); hold on;
pl.LineWidth = 2;
pl = plot(mean(e,2) - std(e,[],2), 'r--'); hold on;
pl = plot(mean(e,2) + std(e,[],2), 'r--'); hold on;


meanE = mean(e,2);
nF = find(meanE == min(meanE));
% nF = 15;

figure;
subplot(1,2,1)
imagesc(nnmfPLV, [0 1]); colorbar;
subplot(1,2,2)
imagesc(W*H, [0 1]); colorbar;


chans = find(~bad);

div = divisors(nF);
nDiv = length(div);
hx = div(nDiv/2); hy = div(nDiv/2 + 1); 


%%
close all

c = distinguishable_colors(nF+1,'k');      

c(1,:) = [0 0 0];


nF = 6;

W = []; H = [];
for it = 1:10
    [W, H] = nnmf(nnmfPLV, nF, 'replicates',10000);
    h = figure('Position', [435 908 1130 413]);


%     for iF = 1:nF
    class = [];
    Wmat = zeros(10,10);
    for ch = 1:length(chans)
        maxF = find(W(ch,:) == max(W(ch,:)));
        Wmat(chMap == chans(ch)) = maxF;
%         Wmat(chMap == chans(ch)) = median(W(ch,iF,it));
        class = [class, maxF];
    end
%     meanW = W(:,:,it);
%         subplot(hx,hy,iF);
%         imagesc(Wmat, [0 max(meanW(:))]);
    nFused = max(unique(Wmat(:)))+1;
    subplot(1,2,1)
    imagesc(Wmat); colormap(c(1:nFused,:)); colorbar;
    axis off
    
    Rg = [];
    for F = 1:nF
%         chF = find(class == F);
%         chNF = find(class ~= F);
%         PLV_F = nnmfPLVo(chF,chF);
%         PLV_NF = nnmfPLVo(chF,chNF);
%         PLV_F(PLV_F==1) = [];
%         PLV_NF(PLV_NF==1) = [];
        [xx, yy] = find(Wmat == F);
        cntr = [mean(xx), mean(yy)]; %center of module
        r = arrayfun(@(X,Y) sqrt((cntr(1) - X)^2) + (cntr(2) - Y)^2, xx, yy) * .4; %distance in mm

        Rg(F) = rms(r);
    end
    subplot(1,2,2)
    plot(1:nF, Rg, 'o'); hold on;
    hline(medan(Rg, 'omitnan'))
%     ylim([    -1 1])
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box off
    ax.YAxis.Label.String = 'PLV within - PLV outside';
    ax.XAxis.Label.String = 'cluster ID';


    
    filename = sprintf('NNMFclustering_iter%02i_wake_v2.pdf', it);
%     savepdf(h, filename)
    
%     end

%      clf;
end
figure;
imagesc(mean(W,3)); colorbar;

[A, I] = sort(distAll);
figure; 
plotPLV = maxPLV(I,I);
% plotPLV(isnan(plotPLV(5,:)), :)= [];
% plotPLV(:, isnan(plotPLV(5,:)))= [];
imagesc(plotPLV, [0 1]); colormap('jet'); colorbar
%%
figure('Position',[1921 127 1000 1000]);
for ch = 1:length(rippleStats.locs)
    chanLabel = str2double(rippleStats.chanLabels{ch});
    
    gridLoc = find(chMap == chanLabel);
    im = zeros(size(chMap));

    for ch2 = 1:length(rippleStats.locs)
        chanLabel2 = str2double(rippleStats.chanLabels{ch2});
        locs2 = rippleStats.locs{ch2};

        gridLoc2 = find(chMap == chanLabel2);

        row = rem(gridLoc2,size(chMap,1));
        if row == 0; row = 10; end
        col = ceil(gridLoc2/size(chMap,1));
        
        
        
        im(col,row) = coOccur(ch,ch2);
%         im(col,row) = maxPLV(ch,ch2);
       
    end
    
    if ~ismember(ch, badChannels)
    
        subplot(size(chMap,1),size(chMap,2),gridLoc);
        imagesc(im,[0 0.6]); hold on;
        colormap('hot')
        axis off
        title(num2str(chanLabel))
    end
%     xlim([-500 500])
    
    
end
%%
% save('T11_PLV_NREM.mat', 'maxPLV', 'PLVstruct', 'coOccur','-v7.3')
save('T11_PLV_NREM.mat', 'maxPLV', 'PLVstruct' ,'-v7.3')



















