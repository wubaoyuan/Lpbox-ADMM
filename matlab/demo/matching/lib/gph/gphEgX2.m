function [Eg1, A1] = gphEgX2(Pt1, X, Eg2, parGph)
% Connecting points to obtain edge.
%
% Input
%   Pt1     -  graph node, 2 x n1
%   X       -  correspondence, n1 x n2
%   Eg2     -  edge connection for the 2nd graph, 2 x m2
%   parGph  -  parameter
%
% Output
%   Eg1     -  edge connection for the 1st graph, 2 x 2m1
%   A1      -  node-node adjacency (has to be symmetric), n0 x n0
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-01-2012

% function parameter
link = ps(parGph, 'link', 'full');
val = ps(parGph, 'val', 3);

% dimension
[n1, n2] = size(X);

% find edge on the original graph
[Eg0, A0] = gphEg(Pt1, parGph);

% find edge
[~, idx2] = asgX2Idx(X);

% change index
ind2 = sub2ind([n1 n1], idx2(Eg2(1, :)), idx2(Eg2(2, :)));
A2 = zeros(n1, n1);
A2(ind2) = 1;
A1 = A0 | A2;

% graph edge
idx = find(triu(A1, 1));
if isempty(idx)
    Eg1 = [];
else
    [is, js] = ind2sub([n1 n1], idx);
    Eg1 = [is, js]';
    Eg1 = [Eg1, Eg1([2 1], :)];
end

%Eg0 = [idx(Eg(1, :)); idx(Eg(2, :))];

