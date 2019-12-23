function [idxTrs, idxTes] = laicSortX(X, Pts, C, H)
% Sort the objective function.
%
% Input
%   X       -  correspondence matrix, Nm x Nt
%   Pts     -  point position, 1 x 2 (cell), 2 x Nm, 2 x Nt
%   H       -  reconstruction matrix of the template point set, Nm x Nm
%   C       -  feature matching cost matrix, Nm x Nt
%
% Output
%   X       -  correspondence matrix, Nm x Nt
%   obj     -  objective
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-13-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-26-2012

% dimension
[Nm, Nt] = size(C);

% objective
[~, ~, ~, obj1s, obj2s] = laicObj(X, Pts, C, H);
objs = obj1s + obj2s;

% sort
[~, idxTrs] = sort(objs);
m = length(idxTrs);
idxTes = zeros(1, m);
for i = 1 : m
    row = idxTrs(i);
    idxTes(i) = find(X(row, :));
end
