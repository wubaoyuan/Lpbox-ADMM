function [Eg0, A0] = gphEgX(Pt0, vis, parGph)
% Connecting points to obtain edge.
%
% Input
%   Pt0     -  graph node, 2 x n0
%   vis     -  outlier index, 1 x n0
%                vis(i) == 0: i is an outlier
%   parGph  -  parameter
%
% Output
%   Eg0     -  edge, 2 x 2m
%   A0      -  node-node adjacency (has to be symmetric), n0 x n0
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-01-2012

% function parameter
link = ps(parGph, 'link', 'full');
val = ps(parGph, 'val', 3);

% dimension
idx = find(vis);
Pt = Pt0(:, vis);
n0 = size(Pt0, 2);
n = size(Pt, 2);

% find edge on the original graph
[Eg0, A0] = gphEg(Pt0, parGph);

% find edge
[Eg, A] = gphEg(Pt, parGph);

% change index
A2 = zeros(n0, n0);
A2(vis, vis) = A;
A0 = A0 | A2;

% graph edge
idx = find(triu(A0, 1));
if isempty(idx)
    Eg = [];
else
    [is, js] = ind2sub([n0 n0], idx);
    Eg = [is, js]';
    Eg = [Eg, Eg([2 1], :)];
end
Eg0 = Eg;

%Eg0 = [idx(Eg(1, :)); idx(Eg(2, :))];

