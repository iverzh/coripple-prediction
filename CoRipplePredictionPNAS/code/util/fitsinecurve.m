function [sFitOrig, modulationOrig, modulation, fvalOrig, fval] = fitsinecurve(xAll,yAll,N)

bins = linspace(-pi,pi,N+1);


y = [];
SEMall = [];
    
for b = 1:length(bins)-1
    ii = xAll > bins(b) & xAll < bins(b+1) & isfinite(yAll);
    y = [y mean(yAll(ii), 'omitnan')];
    SEMall = [SEMall std(yAll(ii), 'omitnan')/sqrt(sum(ii))];
                    
end
x = smoothdata(bins,'movmean',2);
x(1) = [];
    
y = y/max(y);    
yu = max(y);
yl = min(y);
yr = (yu-yl);                               % Range of ?y?
yz = y-yu+(yr/2);
per = 2*pi;                     % Estimate period
ym = mean(y);                               % Estimate offset
fit = @(b,x)  b(1).*(sin(2*pi*x./per + 2*pi/b(2))) + b(3);    % Function to fit
% fit = @(b,x)  b(1).*(cos(2*pi*x./per) + b(3)) + b(2);    % Function to fit
fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
[sFitOrig, fvalOrig] = fminsearch(fcn, [yr; -1; ym]); 
xp = linspace(min(x),max(x));
figure(10)
plot(x,y,'b',  xp,fit(sFitOrig,xp), 'r')
grid

modulationOrig = sFitOrig(1);

fval = [];
modulation = [];
for iter = 1:1000
    yAllShuffle = yAll(randperm(length(yAll)));
    xAllShuffle = xAll(randperm(length(xAll)));
    yShuffle = [];
    for b = 1:length(bins)-1
        ii = xAllShuffle > bins(b) & xAllShuffle < bins(b+1);
        yShuffle = [yShuffle mean(yAllShuffle(ii))];
    
    
    end
%     yShuffle = yShuffle - min(yShuffle);
    
    y = yShuffle/max(yShuffle);
    yu = max(y);
    yl = min(y);
    yr = (yu-yl);                               % Range of ?y?
    yz = y-yu+(yr/2);
    per = 2*pi;                     % Estimate period
    ym = mean(y);                               % Estimate offset
    fit = @(b,x)  b(1).*(sin(2*pi*x./per + 2*pi/-1)) + b(2);    % Function to fit
%     fit = @(b,x)  b(1).*(cos(2*pi*x./per + b(3))) + b(2);    % Function to fit
    fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
    [sFit, fval(iter)] = fminsearch(fcn, [yr;  -1; ym]);                       % Minimise Least-Squares
    xp = linspace(min(x),max(x));
    modulation(iter) = sFit(1);
    
% 
%     if fval(iter) < fvalOrig
%         figure;
%         plot(x,y,'b',  xp,fit(sFit,xp), 'r')
%         grid
%     
%         close;
%     end


    


    
    
end

% figure; 
% histogram(abs(modulation)); hold on
% vline(modulationOrig)


% end

























