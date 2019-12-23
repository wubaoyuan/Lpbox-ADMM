function [KP1, KP, KQ, beta, lamQ] = calKnlGphPQA(gphs, Dess, parKnl)
% Compute node and feature affinity matrix for graph matching.
%
% Remarks
%   The edge feature is asymmetric.
%
% Input
%   gphs    -  graphs, 1 x 2 (cell)
%    G      -  node-edge incidence matrix, ni x mi
%    H      -  node-edge incidence matrix, ni x mi
%   parKnl  -  parameter
%     alg   -  method of computing affinity, {'toy'} | 'pas'
%              'toy': toy data
%              'pas': Pascal data
%              'ucf': shape data
%
% Output
%   KP      -  node-node affinity, n1 x n2
%   KQ      -  edge-edge affinity, m1 x m2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-09-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-16-2012

% function parameter
beta = ps(parKnl, 'beta', 0);
lamQ = ps(parKnl, 'lamQ', 1);
lamP = ps(parKnl, 'lamP', 1);

% dimension
gph1 = gphs{1};
gph2 = gphs{2};
[n1, m1] = size(gph1.G);
[n2, m2] = size(gph2.G);
prIn('calKnlGphPQ', 'n1 %d, n2 %d, m1 %d, m2 %d', n1, n2, m1, m2);

% graph element
[P1, Q1] = stFld(gphs{1}, 'Pt', 'PtD');
[P2, Q2] = stFld(gphs{2}, 'Pt', 'PtD');

% KP (sift)
Dess{1} = double(Dess{1});
Dess{2} = double(Dess{2});
DP = conDstChi(Dess{1}, Dess{2});
KP1 = conKnl(DP, 'nei', .2);

% KP (geometry)
KP2 = 2 * P1' * P2 - (P1 .* P1)' * ones(2, n2) - ones(n1, 2) * (P2 .* P2);
KP2 = KP2 * (1 - lamQ);
KP = lamP * KP1 + KP2;

% KQ
KQ = 2 * Q1' * Q2 - (Q1 .* Q1)' * ones(2, m2) - ones(m1, 2) * (Q2 .* Q2);

% make sure KQ is positive
kqMi = min(KQ(:));
beta = -kqMi + beta;
KQ = KQ + beta;
KQ = KQ * lamQ;
kqMi = min(KQ(:));
kqMa = max(KQ(:));
pr('kqMi %.2f, kqMa %.2f', kqMi, kqMa);

% normalize
% KQ = knlEgNor(KQ, parKnl);

prOut;
