function [Eg, A, vis] = gphEg(Pt, parGph)
% Connecting points to obtain edge.
%
% Input
%   Pt      -  graph node, 2 x n
%   parGph  -  parameter
%     link  -  link type, {'full'} | 'rand' | 'del' | 'del2' | 'knei' | 'enei' | 'non'
%                'full': fully connect
%                'rand': randomly connect
%                'del' : delaunay triangulation
%                'del2': delaunay triangulation with refinement
%                'knei': k nearest neighbor
%                'enei': epsilon nearest neighbor
%                'non' : no connection (no edge)
%     val   -  parameter value, {3}
%                link == 'rand': edge density
%                link == 'knei': value of k (#neighbors)
%                link == 'enei': value of epsilon (the threshold)
%                link == 'del2': minimum #connected nodes for each point
%
% Output
%   Eg      -  edge, 2 x 2m
%   A       -  node-node adjacency (has to be symmetric), n x n
%   vis     -  binary indicator of nodes that have been kept, 1 x n | []
%                only used if link == 'del2'
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 04-26-2012

% function parameter
link = ps(parGph, 'link', 'full');
val = ps(parGph, 'val', 3);

% dimension
n = size(Pt, 2);

% full-connected graph
if strcmp(link, 'full')
    A = ones(n, n);

% randomly connected graph
elseif strcmp(link, 'rand')
    A = rand(n, n) <= val;

% delaunay triangulation
elseif strcmp(link, 'del')
    A = gphAdjDel(Pt);
    
% delaunay triangulation with refinement
elseif strcmp(link, 'del2')
    [A, vis] = gphAdjDelRef(Pt, val);

% k-nearest neighbor
elseif strcmp(link, 'knei')

    % distance
    D = conDst(Pt, Pt);
    D = mdiag(D, inf);

    % neighbor
    [~, Idx] = sort(D);
    A = zeros(n, n);
    for c = 1 : n
        idx1 = sub2ind([n n], zeros(1, val) + c, Idx(1 : val, c)');
        idx2 = sub2ind([n n], Idx(1 : val, c)', zeros(1, val) + c);
        A(idx1) = 1;
        A(idx2) = 1;
    end

% eps-nearest neighbor
elseif strcmp(link, 'enei')
    D = conDst(Pt, Pt);
    D = real(sqrt(D));
    D = mdiag(D, inf);
    A = D <= val;
    
% non edge
elseif strcmp(link, 'non')
    A = [];
    
else
    error('unknown edge link type: %s', link);
end

% make sure A is symmetric and with 0 on diagonal
A = triu(A, 1);
A = A + A';

if ~exist('vis', 'var')
    vis = [];
end

% graph edge
idx = find(triu(A, 1));
if isempty(idx)
    Eg = [];
else
    [is, js] = ind2sub([n n], idx);
    Eg = [is, js]';
    Eg = [Eg, Eg([2 1], :)];
end

% print information
m = length(idx);
prInOut('gphAdj', 'link %s, val %.1f, #points %d, #edges %d', link, val, n, m);
