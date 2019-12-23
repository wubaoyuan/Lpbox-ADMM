function [P, me, sca] = norP(P0)
% Normalize points set.
%
% Input
%   P0      -  original point matrix, d x n
%
% Output
%   P       -  new point matrix, d x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-11-2012

% dimension
[d, n] = size(P0);

% mean
me = mean(P0, 2);

% remove mean
P = P0 - repmat(me, 1, n);

% variance
sca = sqrt(real(max(sum(P .^ 2, 1))));
    
% normalize
P = P / sca;
