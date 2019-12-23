function [KP, KQ] = conKnlGphTra(TraPss, IdxPs, TraQss, IdxQss, P, varargin)
% Compute affinity matrix for trajectory graph matching.
%
% Remark
%   k = k1 x k2
%
% Input
%   TraPss   -  1st order feature, 1 x 2 (cell), 1 x k_i(cell), 2 x n_i^c
%   IdxPs    -  1st order index for available frames, 1 x 2 (cell), 2 x k_i
%   TraQss   -  2nd order feature, 1 x 2 (cell), k_i x k_i (cell), 2 x n_i^c
%   IdxQss   -  2nd order index for available frames, 1 x 2 (cell), 2 x k_i x k_i
%   P        -  time warping, len x 2
%
% Output
%   KP      -  1st order affinity, k1 x k2
%   KQ      -  2nd order affinity (sparse matrix), k x k
%   KE      -  2nd order edge-by-edge affinity, l1 x l2
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 07-10-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
m = length(TraPss);
ks = cellDim(Pss, 2);
k1 = ks(1); k2 = ks(2);
prIn('conKnlGphTra', 'k1 %d, k2 %d', k1, k2);

% 1st order distance
DP = zeros(k1, k2);
for c1 = 1 : k1
    idx1 = IdxPs{1}(:, c1);
    vis1 = P(:, 1) >= idx1(1) & P(:, 1) <= idx1(2);
    
    for c2 = 1 : k2
        idx2 = IdxPs{2}(:, c2);
        vis2 = P(:, 2) >= idx2(1) & P(:, 2) <= idx2(2);
        
        vis = vis1 & vis2;
        pF1s = P(vis, 1) - idx1(1) + 1;
        pF2s = P(vis, 2) - idx2(1) + 1;
        
        X1 = TraPss{1}{c1}(:, pF1s);
        X2 = TraPss{2}{c2}(:, pF2s);
        
        dsts = sum((X2 - X1) .^ 2);
        DP(c1, c2) = sum(sqrt(dsts));
    end
end

% 2nd order
DQ = zeros(k1, k2, k1, k2);
for c1a = 1 : k1
    for c1b = c1a + 1 : k1
        if isempty(TraQss{1}{c1a, c1b})
            continue;
        end
        idx1 = IdxQss{1}(:, c1a, c1b);
        vis1 = P(:, 1) >= idx1(1) & P(:, 1) <= idx1(2);
        
        for c2a = 1 : k2
            for c2b = c2a + 1 : k2
                if isempty(TraQss{2}{c2a, c2b})
                    continue;
                end
                idx2 = IdxQss{2}(:, c2a, c2b);
                vis2 = P(:, 2) >= idx2(1) & P(:, 2) <= idx2(2);
                
                vis = vis1 & vis2;
                pF1s = P(vis, 1) - idx1(1) + 1;
                pF2s = P(vis, 2) - idx2(1) + 1;
               
                X1 = TraQss{1}{c1a, c1b}(:, pF1s);
                X2 = TraQss{2}{c2a, c2b}(:, pF2s);
                dst1s = X1(1, :);
                dst2s = X2(1, :);
                dsts = abs(dst1s - dst2s) ./ (min([dst1s; dst2s], [], 1) + eps);
                    
                ang1s = X1(2, :);
                ang2s = X2(2, :);
                angs = cos(ang1s - ang2s);
                
                vis = dsts < .5 & angs > cos(pi / 8);
                KQ(c1a, c2a, c1b, c2b) = KQ(c1a, c2a, c1b, c2b) + length(find(vis));
            end
        end
    end
end

k = ks(1) * ks(2);
KQ = reshape(KQ, k, k);
KQ = KQ + KQ';

p = sub2ind([k k], 1 : k, 1 : k);
KQ(p) = len;

% store
ali.alg = 'sqm';

prOut;
