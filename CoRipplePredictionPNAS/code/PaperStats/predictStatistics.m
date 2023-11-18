


subjects = {'T11', 'MG29', 'MG63'};
state = 'NREM';

predictRipp = [];
predictNull = [];
sigRipp = [];
sigNull = [];
interact = {};
subject = subjects{1};
stageFolder = sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/SpikePredict/PredictAll/%s/Medial', subject);
load(fullfile(stageFolder, sprintf('SpikeTimes_%s_%s-bspikeRipp_v2-processed-noPreWin.mat', subject, state)))
units = LoadSpikeTimes(subject,'CellExplorer','medial');
SU = find(cellfun(@(X) strcmp(X, 'pyr') | strcmp(X, 'int'), units(:,3)));

for iiA = 1:size(rippStats.predict,1)
    for iiB = 1:size(rippStats.predict,2)
        chA = units{SU(iiA),1}; chB = units{SU(iiB),1};
        typeA = units{SU(iiA),3}; typeB = units{SU(iiB),3};
        if chA ~= chB
            if ~isempty(rippStats.predict{iiA,iiB})
                predictRipp = [predictRipp rippStats.predict{iiA,iiB}];
                predictNull = [predictNull nullStats.predict{iiA,iiB}];

                sigRipp = [sigRipp rippStats.Sig(iiA,iiB)];
                sigNull = [sigNull nullStats.Sig(iiA,iiB)];
            else
                sigRipp = [sigRipp nan];
                predictNull = [predictNull nan];

                predictRipp = [predictRipp nan];
                sigNull = [sigNull 0];
            end
            interact = [interact sprintf('%s-%s', typeA, typeB)];

        end

    end
end

for s = 2:3
    subject = subjects{s};
    stageFolder = sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/SpikePredict/PredictAll/%s', subject);
    load(fullfile(stageFolder, sprintf('SpikeTimes_%s_%s-bspikeRipp_v2-processed-noPreWin.mat', subject, state)))
    units = LoadSpikeTimes(subject,'CellExplorer','medial');
    SU = find(cellfun(@(X) strcmp(X, 'pyr') | strcmp(X, 'int'), units(:,3)));
    
    for iiA = 1:size(rippStats.predict,1)
        for iiB = 1:size(rippStats.predict,2)
            chA = units{SU(iiA),1}; chB = units{SU(iiB),1};
            typeA = units{SU(iiA),3}; typeB = units{SU(iiB),3};
    
            if chA == chA
                if ~isempty(rippStats.predict{iiA,iiB})
                    predictRipp = [predictRipp rippStats.predict{iiA,iiB}];
                    predictNull = [predictNull nullStats.predict{iiA,iiB}];
    
                    sigRipp = [sigRipp rippStats.Sig(iiA,iiB)];
                    sigNull = [sigNull nullStats.Sig(iiA,iiB)];
                else
                    sigRipp = [sigRipp 0];
                    predictNull = [predictNull nan];
    
                    predictRipp = [predictRipp nan];
                    sigNull = [sigNull 0];
                end
                interact = [interact sprintf('%s-%s', typeA, typeB)];

            end
        end
    end
end



%%

ii = cellfun(@(X) strcmp(X, 'pyr-pyr'), interact);
sum(predictRipp > predictNull & (sigRipp == 1 | sigNull == 1) & ii) / sum((sigRipp == 1 | sigNull == 1) & ii)



















