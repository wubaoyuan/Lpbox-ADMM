function gph = newGphCt(Pt, Ct, gphTr)
% Generate a graph by connecting points.
%
% Input
%   Pt      -  graph node, d x n
%   Ct      -  constraints, k x n
%   gphTr   -  training graph
%
% Output
%   gph     -  graph
%     Pt    -  graph node, d x n
%     Eg    -  graph edge, 2 x (2m)
%     G     -  node-edge adjacency, n x m
%     H     -  augumented node-edge adjacency, n x (m + n)
%     PtD   -  displacement, 2 x m
%     dsts  -  distance, 1 x m
%     angs  -  angle, 1 x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-23-2012

% dimension
n = size(Pt, 2);

% edge
Eg = gphEgSame(gphTr.Eg, Ct);

% incidence matrix
[G, H] = gphEg2Inc(Eg, n);

% second-order feature
[PtD, dsts, angs] = gphEg2Feat(Pt, Eg);

% store
gph.Pt = Pt;
gph.Eg = Eg;
gph.G = G;
gph.H = H;
gph.PtD = PtD;
gph.dsts = dsts;
gph.angs = angs;
