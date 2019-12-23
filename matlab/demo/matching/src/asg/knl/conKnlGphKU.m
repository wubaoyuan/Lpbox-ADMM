function K = conKnlGphKU(KP, KQ, gphs)
% Compute the global affinity matrix for graph matching.
%
% Remark
%   nn = n1 x n2
%   The edge is undirected.
%   To deal with directed edge, please use the function "conKnlGphKU.m".
%
% Input
%   KP      -  node-node affinity, n1 x n2
%   KQ      -  edge-edge affinity, m1 x m2
%   gphs    -  graphs, 1 x 2 (cell)
%
% Output
%   K       -  affinity, nn x nn (sparse)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-09-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-07-2013

% symmetric edge feature
KQ = [KQ, KQ; KQ, KQ];

% dimension
[n1, n2] = size(KP);
[m1, m2] = size(KQ);
nn = n1 * n2;
prIn('conKnlGphKU', 'nn %d, m1 %d, m2 %d', nn, m1, m2);

% edge
Eg1 = gphs{1}.Eg;
Eg2 = gphs{2}.Eg;

% global kernel
K = knlPQ2K(KP, KQ, Eg1, Eg2, n1, n2, m1, m2, nn);

prOut;
