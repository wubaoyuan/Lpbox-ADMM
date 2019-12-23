function [KP, KQ] = conKnlGphTra(XPtss, XEgss, Egss, As, P, parKnl)
% Compute affinity matrix for trajectory graph matching.
%
% Remark
%   k0 = k10 x k20
%
% Input
%   XPtss   -  node feature, 1 x 2 (cell), 1 x nFi (cell), dPt_i x k_{t_i}^i
%   XEgss   -  edge feature, 1 x 2 (cell), 1 x nFi (cell), dEg_i x l_{t_i}^i
%   Egss    -  graph edge, 1 x 2 (cell), 1 x nFi (cell), 2 x l_{t_i}^i
%   As      -  trajectory status, 1 x 2 (cell), k0_i x nFi
%   P       -  frame correspondence, nF x 2
%   parKnl  -  parameter, see function conKnlGph for the detailed setting
%     wEg   -  weight for edge feature, [] | dEg x 1
%
% Output
%   KP      -  1st order affinity, k10 x k20
%   KQ      -  2nd order affinity, k0 x k0 (sparse matrix)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 07-10-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function parameter
wEg = ps(parKnl, 'wEg', []);

% dimension
nF = size(P, 1);
k10 = size(As{1}, 1);
k20 = size(As{2}, 1);
k0 = k10 * k20;
prIn('conKnlGphTra', 'k10 %d, k20 %d, nF %d', k10, k20, nF);

% per frame
[KP, CP] = zeross(k10, k20);
KQ = sparse(k0, k0);
CQ = sparse(k0, k0);
prCIn('frame', nF, .1);
for iF = 1 : nF
    prC(iF);
    
    % frame index
    p1 = P(iF, 1);
    p2 = P(iF, 2);

    XPts = {XPtss{1}{p1}, XPtss{2}{p2}};
    XEgs = {XEgss{1}{p1}, XEgss{2}{p2}};
    Egs = {Egss{1}{p1}, Egss{2}{p2}};
    visPts = {As{1}(:, p1), As{2}(:, p2)};
    
    % node & edge affinity
    [KPt, KEg] = knlX2PE(XPts, XEgs, wEg);

    % normalize
    KEg = knlNorEg(KEg, parKnl);

    % convert to graph affinity
    [KPi, KQi] = knlPE2PQ(KPt, KEg, Egs, visPts);
    
    % store
    KP = KP + KPi;
    KQ = KQ + KQi;
    %CP(idxPt1, idxPt2) = CP(idxPt1, idxPt2) + 1;
end
prCOut(nF);

% normalize
%KEg = knlQ2E(KQ);

prOut;
