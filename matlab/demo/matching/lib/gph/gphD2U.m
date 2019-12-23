function gphU = gphD2U(gphD)
% Convert a directed graph to an undirected one.
%
% Input
%   gphD    -  directed graph
%     Pt    -  graph node, d x n
%     Eg    -  graph edge, 2 x mD
%     vis   -  binary indicator of nodes that have been kept, 1 x n | []
%     G     -  node-edge adjacency, n x mD
%     H     -  augumented node-edge adjacency, n x mD
%     PtD   -  edge feature, 2 x mD
%     dsts  -  distance, 1 x mD
%     angs  -  angle, 1 x mD
%
% Output
%   gphU    -  undirected graph
%     Pt    -  graph node, d x n
%     Eg    -  graph edge, 2 x (2m)
%     vis   -  binary indicator of nodes that have been kept, 1 x n | []
%     G     -  node-edge adjacency, n x m
%     H     -  augumented node-edge adjacency, n x (m + n)
%     PtD   -  edge feature, 2 x m
%     dsts  -  distance, 1 x m
%     angs  -  angle, 1 x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% dimensin
n = size(gphD.Pt, 2);

% incidence matrix
[G, H] = gphEg2IncU(gphD.Eg, n);

% store
gphU = gphD;
gphU.G = G;
gphU.H = H;
