function K = conKnlGphKD(KP, KQ, gphs)
% Compute the global affinity matrix for graph matching.
%
% Remarks
%   The edge is directed and the edge feature is asymmetric.   
%   nn = n1 x n2
%
% Input
%   KP      -  node-node affinity, n1 x n2
%   KQ      -  edge-edge affinity, m1 x m2
%   gphs    -  graphs, 1 x 2 (cell)
%
% Output
%   K       -  global affinity, nn x nn (sparse)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-09-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% dimension
[n1, n2] = size(KP);
[m1, m2] = size(KQ);
nn = n1 * n2;
prIn('conKnlGphKD', 'nn %d, m1 %d, m2 %d', nn, m1, m2);

% edge
Eg1 = gphs{1}.Eg;
Eg2 = gphs{2}.Eg;

% global kernel for directed edges
K = knlPQ2K(KP, KQ, Eg1, Eg2, n1, n2, m1, m2, nn);

prOut;
