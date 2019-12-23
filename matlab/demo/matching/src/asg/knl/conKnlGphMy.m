function [Ks, lambs] = conKnlGphMy(KP, KQ, Gs, Hs, egy, parKnl)
% Compute affinity matrix for graph matching.
%
% Remark
%   nnT = n1T x n2T
%   If visPts{1} == [], then visPts{1} is default to ones(ki, 1)
%   If visPts{1} ~= [], then ki = length(find(visPts{i}))
%
% Input
%   XPs       -  node feature, 1 x 2 (cell), dPt x ni
%   XQs       -  graph feature, 1 x 2 (cell), dEg x 2mi
%   Egs       -  graph edges, 1 x 2 (cell), 2 x 2mi
%   parKnl    -  parameter
%     wQ      -  weight for edge feature, [] | dEg x 1
%     visPts  -  node existence status, {[]} | 1 x 2 (cell), niT x 1
%     diag    -  flag of whether setting the diagonal of K as KP, {'y'} | 'n'
%
% Output
%   K         -  affinity, nnT x nnT (sparse)
%   KP        -  node-node affinity, n1T x n2T
%   KQ        -  edge-edge affinity, m1 x m2
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 08-09-2011
%   modify    -  Feng Zhou (zhfe99@gmail.com), 10-16-2011

% function parameter
isDiag = psY(parKnl, 'diag', 'n');

prIn('conKnlGphMy');

% new affinity
G1 = Gs{1};
G2 = Gs{2};
H1 = Hs{1};
H2 = Hs{2};
if isDiag
    L = [KQ, -KQ * G2'; -G1 * KQ, G1 * KQ * G2' + KP];
else
    L = [KQ, -KQ * G2'; -G1 * KQ, G1 * KQ * G2'];    
end

% test affinity
%HH = kron(H2, H1);
%K2 = HH * diag(vec(L)) * HH';
%equal('K', full(K), K2);

% svd
[U, S, V] = svd(L);
s = diag(S);
k = thEgy(s, egy);
pr('ratio: %d/%d %.3f', k, length(s), sum(s(1 : k)) / sum(s));
U = U(:, 1 : k);
V = V(:, 1 : k);
s = s(1 : k);
lambs = s;

% new K
Ks = cell(1, k);
for c = 1 : k
    HH2 = sparse(H2 * diag(V(:, c)) * H2');
    HH1 = sparse(H1 * diag(U(:, c)) * H1');
    Ks{c} = kron(HH2, HH1);

    %Kc2 = s(c) * kron(H2 * diag(V(:, c)) * H2', H1 * diag(U(:, c)) * H1');
end

prOut;
