function [Ps, mes, scas] = norPs(P0s)
% Normalize points set.
%
% Input
%   P0s     -  original point matrix, 1 x m (cell), d x n
%
% Output
%   Ps      -  new point matrix, 1 x m (cell), d x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-16-2012

% dimension
m = length(P0s);

[Ps, mes, scas] = cellss(1, m);
for i = 1 : m
    [Ps{i}, mes{i}, scas{i}] = norP(P0s{i});
end
