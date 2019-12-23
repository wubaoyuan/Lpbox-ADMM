function gphD = gphU2D(gphU)
% Convert an undirected graph to a directed one.
%
% Input
%   gphU    -  undirected graph
%     Pt    -  graph node, d x n
%     Eg    -  graph edge, 2 x mA
%     vis   -  binary indicator of nodes that have been kept, 1 x n | []
%     G     -  node-edge adjacency, n x mA
%     H     -  augumented node-edge adjacency, n x mA
%     PtD   -  edge feature, 2 x mA
%     dsts  -  distance, 1 x mA
%     angs  -  angle, 1 x mA
%
% Output
%   gphD    -  directed graph
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
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-07-2013

% dimensin
n = size(gphU.Pt, 2);

% incidence matrix
[G, H] = gphEg2IncA(gphU.Eg, n);

% store
gphD = gphU;
gphD.G = G;
gphD.H = H;