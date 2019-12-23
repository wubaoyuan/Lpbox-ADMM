function [KP, KQ, beta, lamQ] = conKnlGphPQD(gphs, parKnl)
% Compute node and feature affinity matrix for graph matching.
%
% Remarks
%   The edge is directed and the edge feature is asymmetric.
%
% Input
%   gphs    -  graphs, 1 x 2 (cell)
%    G      -  node-edge incidence matrix, ni x mi
%    H      -  node-edge incidence matrix, ni x mi
%   parKnl  -  parameter
%     alg   -  method of computing affinity, {'toy'} | 'pas' | 'ucf'
%              'toy': toy data
%              'pas': Pascal data
%              'ucf': UCF shape data
%     beta  -  weight to add on the pairwise affinity to make it be possitive, {0}
%     lamQ  -  weight to trade-off between the unary cost and the pairwise cost, {1} 
%
% Output
%   KP      -  node-node affinity, n1 x n2
%   KQ      -  edge-edge affinity, m1 x m2
%   beta    -  weight
%   lamQ    -  weight
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-09-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% function parameter
alg = ps(parKnl, 'alg', 'toy');
beta = ps(parKnl, 'beta', 0);
lamQ = ps(parKnl, 'lamQ', 1);

% dimension
gph1 = gphs{1};
gph2 = gphs{2};
[n1, m1] = size(gph1.G);
[n2, m2] = size(gph2.G);
prIn('conKnlGphPQD', 'alg %s, n1 %d, n2 %d, m1 %d, m2 %d', alg, n1, n2, m1, m2);

% for toy data
if strcmp(alg, 'toy')
    KP = zeros(n1, n2);
    DQ = conDst(gph1.XQ, gph2.XQ);
    KQ = exp(-DQ / .15);

% for Pascal data
elseif strcmp(alg, 'pas')
    DP = conDst(gphs{1}.XP, gphs{2}.XP);
    KP = exp(-real(sqrt(DP)));    

    % distance
    Dst1 = repmat(gph1.dsts', 1, m2);
    Dst2 = repmat(gph2.dsts, m1, 1);
    Dst = abs(Dst1 - Dst2) ./ (min(Dst1, Dst2) + eps);

    % angle
    Ang1 = repmat(gph1.angAs', 1, m2);
    Ang2 = repmat(gph2.angAs, m1, 1);
    Ang = abs(Ang1 - Ang2);
    
    % combine distance and angle
    KQ = exp(-(Dst + Ang) / 2);
    
% for Pascal data
elseif strcmp(alg, 'pas2')
    DP = conDst(gphs{1}.XP, gphs{2}.XP);
    KP = exp(-real(sqrt(DP)));    

    % distance
    Dst1 = repmat(gph1.dsts', 1, m2);
    Dst2 = repmat(gph2.dsts, m1, 1);
    Dst = abs(Dst1 - Dst2) ./ (min(Dst1, Dst2) + eps);

    % angle
    Ang1 = repmat(gph1.angAs', 1, m2);
    Ang2 = repmat(gph2.angAs, m1, 1);
    AngD = Ang1 - Ang2;
    
    AngD2 = abs(AngD);
    
    Vis = AngD < 0;
    AngD(Vis) = -AngD(Vis);
    
    Vis = AngD > pi;
    AngD(Vis) = 2 * pi - AngD(Vis);
    
    % combine distance and angle
    KQ = exp(-(Dst + AngD) / 2);

% for UCF data
elseif strcmp(alg, 'ucf')
    [P1, Q1] = stFld(gphs{1}, 'Pt', 'PtD');
    [P2, Q2] = stFld(gphs{2}, 'Pt', 'PtD');

    % KP
    KP = 2 * P1' * P2 - (P1 .* P1)' * ones(2, n2) - ones(n1, 2) * (P2 .* P2);
    KP = KP * (1 - lamQ);

    % KQ
    KQ = 2 * Q1' * Q2 - (Q1 .* Q1)' * ones(2, m2) - ones(m1, 2) * (Q2 .* Q2);
%    KQ = exp(KQ);

    % make sure KQ is positive
    kqMi = min(KQ(:));
    beta = -kqMi + beta;
    KQ = KQ + beta;
    KQ = KQ * lamQ;
    kqMi = min(KQ(:));
    kqMa = max(KQ(:));
    pr('kqMi %.2f, kqMa %.2f', kqMi, kqMa);
    
% for computational photography Class
elseif strcmp(alg, 'cg')
    XP1 = gphs{1}.XP;
    XP2 = gphs{2}.XP;
    DP = conDstChi(XP1, XP2);
    KP = exp(-real(sqrt(DP)));    

    % distance
    Dst1 = repmat(gph1.dsts', 1, m2);
    Dst2 = repmat(gph2.dsts, m1, 1);
    Dst = abs(Dst1 - Dst2) ./ (min(Dst1, Dst2) + eps);

    % angle
    Ang1 = repmat(gph1.angs', 1, m2);
    Ang2 = repmat(gph2.angs, m1, 1);
    Ang = abs(Ang1 - Ang2);
    
    % combine distance and angle
    KQ = exp(-(Dst + Ang) / 2);
    
else
    error('unknown algorithm: %s', alg);
end

% normalize
% KQ = knlEgNor(KQ, parKnl);

prOut;
