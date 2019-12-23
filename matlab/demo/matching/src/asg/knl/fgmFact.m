function [XQ1, XQ2] = fgmFact(gphs, KP, KQ, beta, lamQ)
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
%
% Output
%   KP      -  node-node affinity, n1 x n2
%   KQ      -  edge-edge affinity, m1 x m2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-09-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-16-2012

% graph element
[P1, Q1, G1, H1] = stFld(gphs{1}, 'Pt', 'PtD', 'G', 'H');
[P2, Q2, G2, H2] = stFld(gphs{2}, 'Pt', 'PtD', 'G', 'H');

% dimension
[n1, m1] = size(G1);
[n2, m2] = size(G2);
prIn('fgmFact', 'n1 %d, n2 %d, m1 %d, m2 %d', n1, n2, m1, m2);

XQ1 = [sqrt(2) * Q1; ones(1, m1); -ones(1, 2) * Q1 .^ 2; ones(1, m1) * sqrt(beta)] * sqrt(lamQ);
XQ2 = [sqrt(2) * Q2; -ones(1, 2) * Q2 .^ 2; ones(1, m2); ones(1, m2) * sqrt(beta)] * sqrt(lamQ);

equal('KQ', KQ, XQ1' * XQ2);

prOut;
