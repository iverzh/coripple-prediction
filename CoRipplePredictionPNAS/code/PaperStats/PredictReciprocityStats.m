

clear
close all
clc

addpath(genpath('../../code'))

%%


subject = 'E1';
state = 'NREM';

units = load(); %see readMe for instructions on how to format unit matrix
matExportFolder = ''; %folder for output of ripple detection.
load(fullfile(matExportFolder,filename));

load('') %load stage mask

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

%% reciprocity B --> A
clc



subjects = {'B1', 'B1-dual', 'E1', 'E2'};
state = 'wake';
clr = [237/255 110/255 87/255; %B1
       255/255 66/255 161/255; %B1-dual
       108/255 227/255 208/255; %E1
       129/255 213/255 83/255]; %E2
% 

BA = []; %B --> A
AB = []; %A --> B

pBA = []; %B --> A
pAB = []; %A --> B

compare = [];
ripp = [];
null = [];
interact = [];
Pr = [];
Pn = [];
countPairs = 0;
figure('Position', [1478 767 318 202]);
for subj =  [1 3] %[1,2,3,4]
    subject = subjects{subj};

    if contains(subject,'dual')
        fdir = sprintf('./out/%s/lateral-medial', 'B1');
        subject = 'B1';
    elseif contains(subject,'B')
        fdir = sprintf('./out/%s/medial', subject);
    else
        fdir = sprintf('./out/%s', subject);
    end
    smoothingWindowExtend = 0;
    win = 150;
    filename = sprintf('Prediction_%s_%s-processed-%i_PreWin-%i_win.mat', subject, state, smoothingWindowExtend, win);
    load(fullfile(fdir, filename))    

    Sr = rippStats.Sig;
    PrTemp = rippStats.predict;
    PrTemp(cellfun(@isempty,PrTemp)) = {nan};
    PrTemp(cellfun(@(X) X > 4, PrTemp)) = {nan};
    PnTemp = nullStats.predict;
    PnTemp(cellfun(@isempty,PnTemp)) = {nan};
    PnTemp(cellfun(@(X) X > 4, PnTemp)) = {nan};
    Z1 = zeros(size(Sr));
    Z2 = zeros(size(Sr));
    n1 = size(Sr,1); % assuming A is the square matrix
    n2 = size(Sr,2); % assuming A is the square matrix
    
    for iA = 1:n1
      for iB = 1:n2
          fOr = rippStats.FilterOut{iA, iB};
          fOn = nullStats.FilterOut{iA, iB};
          countPairs = countPairs + 1;
          [H,p,CI,STATS] = ttest2(fOr, fOn);
          if isempty(H); H =0; CI = 0; end
          if mean(fOr) > mean(fOn) %H && all(CI > 0)
              compare(countPairs) = 1;
          elseif mean(fOr) < mean(fOn) %H && all(CI < 0)
              compare(countPairs) = 2;
          elseif length(CI) > 1
              compare(countPairs) = 3;
          end

          [H,p,CI,STATS] = ttest(fOr);
          if isempty(H); H = 0; CI = 0; end
          ripp(countPairs) = H;

          [H,p,CI,STATS] = ttest(fOn);
          if isempty(H); H = 0; CI = 0; end
          null(countPairs) = H;

          interact{countPairs} = rippStats.CellType{iA, iB};
          Pr(countPairs) = PrTemp{iA,iB};
          Pn(countPairs) = PnTemp{iA,iB};
      end
    end
    
    for iA = 1:n1
      for iB = iA+1:n2
          if ~isnan(PrTemp{iA, iB}) && ~isnan(PrTemp{iB, iA})
              BA = [BA Sr(iA,iB)];
              AB = [AB Sr(iB, iA)];
    
              pBA = [pBA PrTemp{iA,iB}];
              pAB = [pAB PrTemp{iB, iA}];
    
              Z1(iA,iB) = 1;
              Z2(iB,iA) = 1;
          end
    
      end
    end
end
figure;
plot(pAB, pBA, '.')
axis equal

compare(end+1:countPairs) = 0;
meetCrit = sum(compare > 0);
fprintf('%.4f (%i/%i) cell pairs meet inclusion criteria\n', meetCrit/countPairs, meetCrit, countPairs)
fprintf('%.4f (%i/%i) cell pairs are greater during ripples\n', sum(compare == 1)/meetCrit, sum(compare == 1), meetCrit)
fprintf('%.4f (%i/%i) cell pairs are greater during null\n', sum(compare == 2)/meetCrit, sum(compare == 2), meetCrit)
fprintf('%.4f (%i/%i) cell pairs are greater during neither\n', sum(compare == 3)/meetCrit, sum(compare == 3), meetCrit)

s = 'pyr->pyr';
iiType = strcmp(interact, s); % & (null | ripp);
mu = mean(Pr(iiType) - Pn(iiType), 'omitnan');
[H,p,CI,STATS] = ttest(Pr(iiType) - Pn(iiType));
fprintf('%s coR - noR difference: %.4f (%.4f - %.4f) p = %.10f\n', s, mu, CI(1), CI(2), p)
fprintf('%.4f (%i/%i) %s cell pairs are greater during ripples\n', sum(compare(iiType) == 1)/sum(iiType), sum(compare(iiType) == 1), sum(iiType), s)
fprintf('%.4f (%i/%i) %s cell pairs are greater during null   \n', sum(compare(iiType) == 2)/sum(iiType), sum(compare(iiType) == 2), sum(iiType), s)
s = 'pyr->int';
iiType = strcmp(interact, s); % & (null | ripp);
mu = mean(Pr(iiType) - Pn(iiType), 'omitnan');
[H,p,CI,STATS] = ttest(Pr(iiType) - Pn(iiType));
fprintf('%s coR - noR difference: %.4f (%.4f - %.4f) p = %.10f\n', s, mu, CI(1), CI(2), p)
fprintf('%.4f (%i/%i) %s cell pairs are greater during ripples\n', sum(compare(iiType) == 1)/sum(iiType), sum(compare(iiType) == 1), sum(iiType), s)
fprintf('%.4f (%i/%i) %s cell pairs are greater during null   \n', sum(compare(iiType) == 2)/sum(iiType), sum(compare(iiType) == 2), sum(iiType), s)

s = 'int->pyr';
iiType = strcmp(interact, s); % & (null | ripp);
mu = mean(Pr(iiType) - Pn(iiType), 'omitnan');
[H,p,CI,STATS] = ttest(Pr(iiType) - Pn(iiType));
fprintf('%s coR - noR difference: %.4f (%.4f - %.4f) p = %.10f\n', s, mu, CI(1), CI(2), p)
fprintf('%.4f (%i/%i) %s cell pairs are greater during ripples\n', sum(compare(iiType) == 1)/sum(iiType), sum(compare(iiType) == 1), sum(iiType), s)
fprintf('%.4f (%i/%i) %s cell pairs are greater during null   \n', sum(compare(iiType) == 2)/sum(iiType), sum(compare(iiType) == 2), sum(iiType), s)

s = 'int->int';
iiType = strcmp(interact, s); % & (null | ripp);
mu = mean(Pr(iiType) - Pn(iiType), 'omitnan');
[H,p,CI,STATS] = ttest(Pr(iiType) - Pn(iiType));
fprintf('%s coR - noR difference: %.4f (%.4f - %.4f) p = %.10f\n', s, mu, CI(1), CI(2), p)
fprintf('%.4f (%i/%i) %s cell pairs are greater during ripples\n', sum(compare(iiType) == 1)/sum(iiType), sum(compare(iiType) == 1), sum(iiType), s)
fprintf('%.4f (%i/%i) %s cell pairs are greater during null   \n', sum(compare(iiType) == 2)/sum(iiType), sum(compare(iiType) == 2), sum(iiType), s)
% Compute the correlation coefficient and p-value
[r, p] = corr(pAB', pBA', 'type', 'pearson', 'rows','pairwise');

% Display the results
disp(['Pearson correlation coefficient = ' num2str(r)]);
disp(['P-value = ' num2str(p)]);


% Define the data
observed = [sum(AB == 0 & BA == 0) sum(AB == 1 & BA == 0);
            sum(AB == 0 & BA == 1) sum(AB == 1 & BA == 1)]; % observed frequencies

% Compute the expected frequencies
rowTotals = sum(observed, 2); % row totals
colTotals = sum(observed, 1); % column totals
n1 = sum(rowTotals);
expected = rowTotals * colTotals / n1;

% Compute the chi-squared statistic
chiSquared = sum(sum((observed - expected).^2 ./ expected));

% Compute the degrees of freedom and p-value
df = (size(observed, 1) - 1) * (size(observed, 2) - 1);
pValue = 1 - chi2cdf(chiSquared, df);

% Display the results
disp(['Chi-squared statistic = ' num2str(chiSquared)]);
disp(['Degrees of freedom = ' num2str(df)]);
disp(['P-value = ' num2str(pValue)]);
%% reciprocity A --> B B --> C C --> A
AyByCyes = []; %B --> A
AyByCno = []; %A --> B
AyByCyesOther = []; %B --> A
AyByCnoOther = []; %A --> B
AyBnCyes = []; %B --> A
AyBnCno = []; %A --> B
AnByCyes = []; %B --> A
AnByCno = []; %A --> B
AnBnCyes = []; %B --> A
AnBnCno = []; %A --> B
countTriplets =  1;

clc
for subj =  [1,3] %, 4]% [1, 3,4]
    subject = subjects{subj};

    if contains(subject,'T')
        fdir = sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/SpikePredict/PredictAll/%s/Medial', subject);
    else
        fdir = sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/SpikePredict/PredictAll/%s', subject);
    end
    load(fullfile(fdir, sprintf('SpikeTimes_%s_%s-bspikeRipp_v2-processed-noPreWin.mat', subject, state)))
    
    Sr = rippStats.Sig;
    Sn = nullStats.Sig;
    PrTemp = rippStats.predict;
    PrTemp(cellfun(@isempty,PrTemp)) = {nan};
    PrTemp(cellfun(@(X) X > 2, PrTemp)) = {nan};
    Z1 = zeros(size(Sr));
    Z2 = zeros(size(Sr));
    n1 = size(Sr,1); % assuming A is the square matrix
    
    ABC = [];
    null = [];
    for iA = 1:n1
      for iB = 1:n1
          for iC = 1:n1
              for iD = 1
            
    
                  if ~isnan(PrTemp{iA, iB}) && ~isnan(PrTemp{iB, iC}) &&  ~isnan(PrTemp{iC, iA}) && length(unique([iA iB iC])) == 3
                      
%                       ABC = [ABC sum([P{iA, iB} P{iB, iC} P{iC, iA}])];
%                       null = [null ,P{iA, iB} ,P{iB, iC} ,P{iC, iA}];
        
                      
                      if     logical(Sr(iA, iB)) && logical(Sr(iB, iC)) && logical(Sr(iC, iA)) && ~ismember(iC, [iA iB])%B --> C
                          AyByCyes = [AyByCyes 1];
                      elseif logical(Sr(iA, iB)) && logical(Sr(iB, iC)) && ~logical(Sr(iC, iA)) && ~ismember(iC, [iA iB])%B --> C
                          AyByCno  = [AyByCno 1];
                      end

                      if     logical(Sr(iA, iB)) && logical(Sr(iB, iC)) && ~ismember(iC, [iA iB])%B --> C
                          C = Sr(iC,:);
                          C([iA iB]) = [];
                          AyByCyesOther = [AyByCyesOther sum(C)];
                          AyByCnoOther  = [AyByCnoOther sum(C == 0)];
                      end
            
                      if     logical(Sr(iA, iB)) && ~logical(Sr(iB, iC)) && logical(Sr(iC, iA)) && ~ismember(iC, [iA iB])%B --> C
                          AyBnCyes = [AyBnCyes 1];
                      elseif logical(Sr(iA, iB)) && ~logical(Sr(iB, iC)) && ~logical(Sr(iC, iA)) && ~ismember(iC, [iA iB])%B --> C
                          AyBnCno  = [AyBnCno 1];
                      end
        
                      if     ~logical(Sr(iA, iB)) && logical(Sr(iB, iC)) && logical(Sr(iC, iA)) && ~ismember(iC, [iA iB])%B --> C
                          AnByCyes = [AnByCyes 1];
                      elseif ~logical(Sr(iA, iB)) && logical(Sr(iB, iC)) && ~logical(Sr(iC, iA)) && ~ismember(iC, [iA iB])%B --> C
                          AnByCno  = [AnByCno 1];
                      end
        
                      if     ~logical(Sr(iA, iB)) && ~logical(Sr(iB, iC)) && logical(Sr(iC, iA)) && ~ismember(iC, [iA iB])%B --> C
                          AnBnCyes = [AnBnCyes 1];
                      elseif ~logical(Sr(iA, iB)) && ~logical(Sr(iB, iC)) && ~logical(Sr(iC, iA)) && ~ismember(iC, [iA iB])%B --> C
                          AnBnCno  = [AnBnCno 1];
                      end
    
                      countTriplets =  countTriplets + 1;
                  end
              end
              
            
    
                  
              
          end
    
      end
    end
end
clc
fprintf('A--> B yes  B--> C yes  C --> A %.4f || %.4f of all\n', length(AyByCyes) / (length(AyByCyes) + length(AyByCno)), length(AyByCyes) / countTriplets)
fprintf('A--> B yes  B--> C no   C --> A %.4f || %.4f of all\n', length(AyBnCyes) / (length(AyBnCyes) + length(AyBnCno)),  length(AyBnCyes) / countTriplets)
fprintf('A--> B no   B--> C yes  C --> A %.4f || %.4f of all\n', length(AnByCyes) / (length(AnByCyes) + length(AnByCno)),  length(AnByCyes) / countTriplets)
fprintf('A--> B no   B--> C no   C --> A %.4f || %.4f of all\n', length(AnBnCyes) / (length(AnBnCyes) + length(AnBnCno)),  length(AnBnCyes) / countTriplets)


% observed = [sum(AyByCyes) sum(AyBnCyes) sum(AnByCyes) sum(AnBnCyes);
%             sum(AyByCno)  sum(AyBnCno)  sum(AnByCno)  sum(AnBnCno)]; % observed frequencies
observed = [sum(AyByCyes) sum(AyByCyesOther);
            sum(AyByCno)  sum(AyByCnoOther)]; % observed frequencies

% Compute the expected frequencies
rowTotals = sum(observed, 2); % row totals
colTotals = sum(observed, 1); % column totals
n1 = sum(rowTotals);
expected = rowTotals * colTotals / n1;

% Compute the chi-squared statistic
chiSquared = sum(sum((observed - expected).^2 ./ expected));

% Compute the degrees of freedom and p-value
df = (size(observed, 1) - 1) * (size(observed, 2) - 1);
pValue = 1 - chi2cdf(chiSquared, df);

% Display the results
disp(['Chi-squared statistic = ' num2str(chiSquared)]);
disp(['Degrees of freedom = ' num2str(df)]);
disp(['P-value = ' num2str(pValue)]);
% ABCshuffle = [];
% for ii = 1:length(ABC)
%     Y = randsample(null, 3);
%     ABCshuffle = [ABCshuffle, sum(Y)];
% end
figure('Position', [992 904 154 417]); 
bar(1, length(AyByCyes) / (length(AyByCyes) + length(AyByCno))); hold on;
bar(2, length(AnBnCyes) / (length(AnBnCyes) + length(AnBnCno))); hold on;

figure('Position', [992 904 154 417]); 
bar(1, sum(AyByCyes) / (sum(AyByCyes) + sum(AyByCno))); hold on;
bar(2, sum(AyByCyesOther) / (sum(AyByCyesOther) + sum(AyByCnoOther))); hold on;
