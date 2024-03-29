function [b b0] = circ_skewness(alpha, w)

% [b b0] = circ_skew(alpha,w)
%   Calculates a measure of angular skewness.
%
%   Input:
%     alpha     sample of angles
%     [w        weightings in case of binned angle data]
%
%   Output:
%     b         skewness (from Pewsey)
%     b0        alternative skewness measure (from Fisher)
%
%   References:
%     Pewsey, Metrika, 2004
%     Statistical analysis of circular data, Fisher, p. 34
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de

alpha = alpha(:);

if nargin < 2 
  w = ones(size(alpha));
else 
  w = w(:);
end

% compute neccessary values
R = circ_r(alpha,w);
theta = circ_mean(alpha,w);
[foo rho2 mu2] = circ_moment(alpha,w,2,true);

% compute skewness 
b = w'*(sin(2*(circ_dist(alpha,theta))))/sum(w);
b0 = rho2*sin(circ_dist(mu2,2*theta))/(1-R)^(3/2);    % (formula 2.29)


