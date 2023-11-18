
close all
clear 

addpath(genpath('./ASR-master'))
%% Inputs

removalMethod = 'template'; % 'ASR', 'template', 'clip'
fsRipple = 1000; % desired downsample rate
fs = 30e3; %sampling rate
startTime = 15; %samples
endTime = 50; %samples
recenterWin = 60;

% rawData is the raw high sampling rate data. nChannels x samples
rawData = [];

% units is a cell array that is nUnits x 2
% the first column contains the channel number that the unit is detected on
% the second column contains an array with the spike times in seconds
units = {};

%% Remove unit spikes from LFP
dataUnitsRemoved = rawData;
% clear rawData %may be necesarry for memory

for unit = 1:size(units,1)
    
 
    chanNum = units{unit,1};

    dataChan = double(dataUnitsRemoved(chanNum,:));
    spikeTimes = units{unit,2} * fs;
    spikeTimes(spikeTimes > (length(dataChan) - recenterWin)) = [];
    spikeTimes(spikeTimes < recenterWin) = [];
    nSpk  = length(spikeTimes);
    
    if nSpk > 0


        if strcmp(removalMethod,'ASR')
            
            hp.Fs       = fs; %sampling rate Hz
            hp.a1ms     = 225; %time (msec.) of extraction after the spike trough
            hp.b1ms     = 225; %time (msec.) of extraction before the spike trough
            hp.Fpl      = 2;%lower limit for pspectrum to find the minimum frequency of decomposition
            hp.Fpu      = 200;%upper limit for pspectrum to find the minimum frequency of decomposition
            hp.MAIa     = 5; %miss-alignment time interval (msec.) after the current spike indices
            hp.MAIb     = 5; %miss-alignment time interval (msec.) before the current spike indices 
            hp.l_cutoff = 2; %setting a cut off frequency for removing low freq. decomposition depending on how noisy the data is

            in = [];
            spikeTimes(spikeTimes < (hp.b1ms * (fs/1e3) + fs/20)) = [];
            spikeTimes(spikeTimes > (length(dataChan) - (hp.a1ms * (fs/1e3)))) = [];

            in.Spike_cell{1} = round(spikeTimes);
            in.Spike_trial = 1;
            in.data_wb{1} = dataChan;

            [hp,in] = Spike_align(hp,in);

            out = ASR_Func(hp,in); 
            
            dataChanRemUnits = out.data_wb{1};
            
            

        else
        
            spikeMat = zeros(nSpk,startTime+endTime+1);
            for spk = 1:nSpk
                times = round(spikeTimes(spk)-recenterWin:spikeTimes(spk)+recenterWin); %times for recentering 

                spikeOrig = dataChan(times);
                shift = find(spikeOrig == min(spikeOrig), 1, 'first'); 
                shift = shift - round((2*recenterWin+1)/2); %make sure spike times are aligned at peak
            %     
                spikeTimes(spk) = spikeTimes(spk) + shift;

                times = round(spikeTimes(spk)-startTime:spikeTimes(spk)+endTime);
                spikeRecenter = dataChan(times);
                spike = spikeRecenter - spikeRecenter(1);


                switch removalMethod
                    case 'clip'
                        spikeCut = linspace(spikeRecenter(1),spikeRecenter(end),length(spikeRecenter));
                        spike = spikeCut;
                        dataChan(times) = spike;
                end

                spikeMat(spk,:) = spike;
            end

            remSpikes = mean(spikeMat);
            stdSpikes = std(spikeMat);
        
            dataChanRemUnits = dataChan;

            switch removalMethod 
                case 'template'
                    for spk = 1:nSpk
                        times = round(spikeTimes(spk)-startTime:spikeTimes(spk)+endTime);

                        spikeOrig = dataChan(times);
                        spike = spikeOrig;
                        spike = spike  - remSpikes; %subtract template
                        spikeMat(spk,:) = spike;
                        dataChanRemUnits(times) = spike;
                    end
            end

        end

        dataUnitsRemoved(chanNum,:) = dataChanRemUnits;
        fprintf('unit %i cut out (channel %i) \n', unit, chanNum)
    end
end

%% Downsample 

LFP = [];
count = 1;
for ch = 1:size(dataUnitsRemoved,1)
    tmp = double(dataUnitsRemoved(ch,:));

    
    LFP(count,:) = resample(tmp,fsRipple,fs);
    
    % notch filter data to remove line noise (60 Hz + harmonics)
    for nf = 60:60:fsRipple/2 %240
        Wo = nf/(fsRipple/2);
        BW = Wo/35;
        [b,a] = iirnotch(Wo, BW);
        LFP(count,:) = filtfilt(b,a,LFP(count,:));
    end

    fprintf([num2str(ch),' done \n'])
    count = count + 1;
    
end