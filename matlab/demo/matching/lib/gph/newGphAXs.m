function gphs = newGphAXs(Pts, X, parGph)
% Connect nodes to generate edges.
%
% Input
%   Pts     -  graph node, 1 x mG (cell), d x ni
%   X       -  correspondence, n1 x n2
%   parGph  -  parameter for computing the adjacency matrix
%              see gphEg.m for more details
%
% Output
%   gphs    -  graph, 1 x mG (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-01-2012

% dimension
mG = length(Pts);
[n1, n2] = size(X);

% index
[idx1, idx2] = asgX2Idx(X);
idxs = {idx1, idx2};

% per graph
gphs = cell(1, mG);
for iG = 1 : mG
    vis = idxs{iG} ~= 0;
    gphs{iG} = newGphAX(Pts{iG}, vis, parGph);
end
