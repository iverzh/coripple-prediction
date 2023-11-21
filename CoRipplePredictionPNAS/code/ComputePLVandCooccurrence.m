
clear
close all
clc
addpath(genpath('../../code'))

% compute co-ripple rate and phase locking within utah arrays
%%



subject = 'E1';
minOverlap = 25; %ms
plvSmoothWin = 10; %ms
win =  500; %ms 
randWin = [-10000 -2000];
nShuff = 100;
f = ''; %path to LFP data
load(f, 'data', 'chan_labels');
LFP = data; clear data;

matExportFolder = ''; %folder for output of ripple detection.
load(fullfile(matExportFolder,filename));

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
PLVstruct = [];
nCh = size(RBphaseAll,1);
maxPLV = nan(nCh,nCh);
coOccur = nan(nCh,nCh);

for chA = 1:nCh
    for chB = 1:nCh
        clc
        fprintf('calculating plv for %i , %i \n', chA, chB)
        
        rippOverlap = rippMask(chA,:) & rippMask(chB,:);
        b = mask2bounds(rippOverlap);
        overlapDuration = b(:,2) - b(:,1);
        b = b(overlapDuration > minOverlap,:);
        coOccurRate = size(b,1) / length(rippleStats.locs{chA});

        coOccur(chA,chB) = coOccurRate;

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


















