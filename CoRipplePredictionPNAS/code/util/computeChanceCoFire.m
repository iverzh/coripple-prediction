
function coFireAll = computeChanceCoFire(xa, xb, mask, coFireThresh, sO, nIter)


stageInd = find(mask);
xa = find(ismember(stageInd, round(xa)));
xb = find(ismember(stageInd, round(xb)));


% ISIb = diff(xb);
ISIb = [ 0 diff(xb)];

ISIa = diff(xa);
coFireAll = [];
% figure;
for n = 1:nIter
    shuffi  = randperm(length(ISIb),length(ISIb));
    ISIshuff = ISIb(shuffi);
%     xbShuff = randperm(length(stageInd),length(xb));
    xbShuff = xb(1);
    for b = 1:length(ISIshuff); xbShuff = [xbShuff xbShuff(end) + ISIshuff(b)]; end
    
%     shuffi  = randperm(length(ISIa),length(ISIa));
%     ISIshuff = ISIa(shuffi);
% %     xbShuff = randperm(length(stageInd),length(xb));
    xaShuff = xa;
%     for a = 1:length(ISIshuff); xaShuff = [xaShuff xaShuff(end) + ISIshuff(a)]; end


    s = xaShuff - xbShuff';
% 
%     No = histcounts(sO, -50.5:50.5);
%     plot(-50:50, smoothdata(No, 'gaussian', 10)); hold on
    N = histcounts(s, -50.5:50.5);
%     plot(-50:50, smoothdata(N, 'gaussian', 10)); hold on
%     vline(-11:11:11)
%     ylim([0 max(No)])
% 
%     waitforbuttonpress; clf;

%     coF = sum(s(:) <= coFireThresh & s(:) >= 0);

    coFireAll(n,:) = N;
end
% 
% b = mask2bounds(mask);
% dur  = (b(:,2) - b(:,1));
% XX = []; for ss = 1:length(dur); XX(ss) = sum(dur(1:ss)); end
% XX(2:end+1) = XX;
% XX(1) = 0;
% coFireAll = [];
% for n = 1:nIter
%     coF = 0;
%     for ii = 2:length(XX)
%         xaEvent = xa(xa > XX(ii-1) & xa <= XX(ii));
%         xbEvent = xb(xb > XX(ii-1) & xb <= XX(ii));
% 
%         if ~isempty(xaEvent) && ~isempty(xbEvent)
%             xaShuffle = randperm((XX(ii)-XX(ii-1)), length(xaEvent));
%             xbShuffle = randperm((XX(ii)-XX(ii-1)), length(xbEvent));
% 
%             s = xaShuffle - xbShuffle';
%             s = abs(s);
% 
%             coF = coF + sum(s(:) <= coFireThresh & s(:) >= 0);
% 
%         end
%     end
% 
% 
% 
%     coFireAll(n) = coF;
% end

return



%%
% close all
% 
% 
figure; 
yyaxis right
plot(xa, ones([1,length(xa)]), 'b.'); hold on;
plot(xb, 2*ones([1,length(xb)]), 'r.'); hold on;
ax = gca;
% ax.XTick = 0:5:ax.XTick(end)
ylim([0 3])

yyaxis left
plot(sum(rippMask))
% vline(XX)
%%
figure; 
yyaxis right
plot(uTcoRa, ones([1,length(uTcoRa)]), 'b.'); hold on;
plot(uTcoRb, 2*ones([1,length(uTcoRb)]), 'r.'); hold on;
ax = gca;
% ax.XTick = 0:5:ax.XTick(end)
ylim([0 3])

yyaxis left
% plot(sum(rippMask))
plot(coRipMask, 'r-'); hold on;
plot(0.8*coRipMaskO, 'b-'); hold on;
ylim([0 2])
%%

uTb = units{uB,2} * 1e3;
uTb1 = units{uB+1,2} * 1e3;
uTb2 = units{uB+2,2} * 1e3;
figure; 
yyaxis right
plot(uTa, ones([1,length(uTa)]), 'b.'); hold on;
plot(uTb, 2*ones([1,length(uTb)]), 'r.'); hold on;
% plot(uTb1, 3*ones([1,length(uTb1)]), 'r.'); hold on;
% plot(uTb2, 4*ones([1,length(uTb2)]), 'r.'); hold on;
ax = gca;
% ax.XTick = 0:5:ax.XTick(end)
ylim([0 3])

yyaxis left
% plot(sum(rippMask))
plot(coRipMask, 'r-'); hold on;
plot(0.8*coRipMaskO, 'b-'); hold on;
ylim([0 2])

% xlim([5.220e7 5.222e7]) %T11 wake

%% 
% c = 0;
% for iX = 2:length(XX)
%     aT = xa > XX(iX-1) & xa <= XX(iX);
%     bT = xb > XX(iX-1) & xb <= XX(iX);
% 
%     if sum(aT) > 0 && sum(bT) > 0
%         c = c + 1;
%         figure; 
%         
%         plot(xa(aT), ones([1,sum(aT)]), 'b.'); hold on;
%         plot(xb(bT), 2*ones([1,sum(bT)]), 'r.'); hold on;
%         
%         ax = gca;
%         ax.XTick = 0:5:ax.XTick(end)
%         ylim([0 3])
%         grid on;
%         pause; close;
% 
%     end
%     
% 



end










