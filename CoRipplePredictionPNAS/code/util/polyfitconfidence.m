function [p, S, CI ]= polyfitconfidence(x,y)


fits = [];
for iter = 1:1000
    Y = datasample(1:length(x),length(x));
    [pShuff, Sshuff] = polyfit(x(Y)',y(Y)',1);    
    
    x1 = linspace(0,max(x));
    [y1,delta] = polyval(pShuff,x1,Sshuff);

    fits(iter,:) = y1;



end

CI = quantile(fits, [0.025, 0.975]);


[p, S] = polyfit(x',y',1);
% [y1,delta] = polyval(p,x1,S);


return