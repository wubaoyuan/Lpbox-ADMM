function [KU, KQU] = knlGphKD2U(KP, KQ, gphUs)
% Compute the global affinity matrix for undirected graph matching.
%
% Remark
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

m1a = round(m1 / 2);
m2a = round(m2 / 2);

% KQU = KQ(1 : m1a, 1 : m2a);
% KQU2 = [KQU, KQU; KQU, KQU];

KQU = (KQ(1 : m1a, 1 : m2a) + KQ(1 : m1a, m2a + 1 : end)) / 2;
KQU2 = [KQU, KQU; KQU, KQU];

% edge
Eg1 = gphUs{1}.Eg;
Eg2 = gphUs{2}.Eg;

% global kernel for asymmetric edges
KU = knlPQ2K(KP, KQU2, Eg1, Eg2, n1, n2, m1, m2, nn);
