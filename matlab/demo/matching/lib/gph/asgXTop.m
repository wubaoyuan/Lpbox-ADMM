function [X, Y] = asgXTop(X0, KP, m)
% Find the best matching.
%
% Input
%   X0      -  original correspondence matrix, n1 x n2
%   KP      -  node affinity matrix, n1 x n2
%   m       -  #best pairs
%
% Output
%   X       -  new correspondence matrix, n1 x n2
%   Y       -  rank matrix, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-12-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 12-05-2012

% dimension
[n1, n2] = size(X0);

% original pairs
idx = find(X0(:));
m0 = length(idx);

% check the value
vals = KP(idx);
[~, rank] = sort(vals, 'descend');

Y = zeros(n1, n2);
Y(idx(rank)) = 1 : m0;
X = double(Y ~= 0 & Y <= m);