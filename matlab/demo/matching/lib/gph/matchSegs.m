function [segs, T0s, Ts, Ps] = matchSegs(Xs, seg0s, paraM, varargin)
% Find the optimum matching across segmentations of several sequences.
%
% Global
%   K       -  kernel matrix, n x n
%
% Input
%   Xs      -  set of sample matrices, 1 x nS (cell), dim x nFi
%   seg0s   -  set of segmentations, 1 x nS (cell)
%   paraM
%     alg   -  algorithm for finding the matching, 'hug' | 'bb'
%     kM    -  #clusters
%     sim   -  similarity measure between segments, 'dtak' | 'spat' | 'over' | 'hist'
%              'dtak': Dynamic Time Alignment Kernel (DTAK)
%              'spat': spatial distance
%              'hist': Histogram
%              'over': Overlapping of different segmentations on the same sequence. In this case, Xs is useless.
%   varargin
%     save option
%
% Output
%   segs    -  segmentation after matching, 1 x nSub (cell)
%   T0s     -  confusion matrix before adjusting, k x k x nSub x nSub
%   Ts      -  confusion matrix after adjusting, k x k x nSub x nSub
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-18-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
[alg, k, sim] = stFld(paraM, 'alg', 'kM', 'sim');

% save option
[svL, path] = psSv(varargin, 'subx', ['match_' alg '_' sim]);

% load
if svL == 2 && exist(path, 'file')
    prom('m', 'old match segs: %s %s\n', alg, sim);
    [segs, T0s, Ts, Ps] = matFld(path, 'segs', 'T0s', 'Ts', 'Ps');
    return;
end
prom('m', 'new match segs: %s %s...', alg, sim);

% insert empty clusters
nS = length(Xs);
for iS = 1 : nS
    Gi0 = seg0s{iS}.G;
    [ki, mi] = size(Gi0);
    if ki < k
        seg0s{iS}.G = [Gi0; zeros(k - ki, mi)];
    end
end

% segment similarity matrix
if strcmp(sim, 'dtak')
    T0s = genConfDtak(Xs, seg0s);
elseif strcmp(sim, 'spat')
    T0s = genConfSpace(Xs, seg0s);
elseif strcmp(sim, 'over')
    T0s = genConfOver(seg0s);
elseif strcmp(sim, 'hist')
    T0s = genConfHist(Xs, seg0s);
else
    error(['unknown similarity type: ' sim]);
end

% multi-dimensional assignment
[Ps, Ts] = asgM(T0s, alg);

% adjust segmentation
segs = cell(1, nS);
segs{1} = seg0s{1};
for i = 2 : nS
    P = Ps(:, :, 1, i);
    G0 = seg0s{i}.G;

    idx = zeros(1, k);
    for c = 1 : k
        p = P(c, :);
        idx(c) = find(p);
    end
    G = G0(idx, :);

    segs{i} = seg0s{i};
    segs{i}.G = G;
end

prom('m', '\n');

% save
if svL > 0
    save(path, 'segs', 'T0s', 'Ts', 'Ps');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ts = genConfDtak(Xs, segs)
% Compute the class confusion matrix for pair of segmentation based on the temporal similarity.
%
% Input
%   Xs    -  set of sample matrix, 1 x nS (cell), dim x nF
%   segs  -  set of segmentation, 1 x nS (cell)
%
% Output
%   Ts    -  confusion matrix, k x k x nS x nS

global KG;

nS = length(Xs);
k = size(segs{1}.G, 1);

% frame lengths
nFs = cellDim(Xs, 2);
s = n2s(nFs);

% frame similarity
% X = cat(2, Xs{:});
% K = conKnl(conDist(X, X), 'nei', nei);

Ts = zeros(k, k, nS, nS);
for iS = 1 : nS
    for jS = iS : nS

        % segment similarity
        Kij = KG(s(iS) : s(iS + 1) - 1, s(jS) : s(jS + 1) - 1);
        segi = segs{iS};
        segj = segs{jS};
        sS = dtaksFord(Kij, segi, segj);

        % sum of pairwise similarity
        Gi = segi.G;
        Gj = segj.G;
        Ts(:, :, iS, jS) = Gi * sS * Gj'; 
        Ts(:, :, jS, iS) = Ts(:, :, iS, jS)';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cs = genConfSpace(Xs, segs)
% Compute the class confusion matrix for pair of segmentation based on the subspace angles.
%
% Input
%   Xs    -  set of sample matrix, 1 x nS (cell), dim x nF
%   segs  -  set of segmentation, 1 x nS (cell)
%
% Output
%   Cs    -  confusion matrix, k x k x nS x nS

nS = length(Xs);
k = size(segs{1}.G, 1);

Cs = zeros(k, k, nS, nS);
for iS = 1 : nS
    for jS = iS + 1 : nS
        Xi = Xs{iS}; segi = segs{iS}; Ii = seg2IH(segi);
        Xj = Xs{jS}; segj = segs{jS}; Ij = seg2IH(segj);
        
        for ci = 1 : k
            for cj = 1 : k
                Yi = Xi(:, Ii(ci, :) == 1);
                Yj = Xj(:, Ij(cj, :) == 1);
        
                D = conDistSpace({Yi, Yj}, 'alg', 'ang');
                Cs(ci, cj, iS, jS) = D(1, 2);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cs = genConfHist(Xs, segs)
% Compute the class confusion matrix for pair of segmentation based on the histogram distance.
%
% Input
%   Xs    -  set of sample matrix, 1 x nS (cell), dim x nF
%   segs  -  set of segmentation, 1 x nS (cell)
%
% Output
%   Cs    -  confusion matrix, k x k x nS x nS

nS = length(Xs);
k = size(segs{1}.H, 1);

% construct the histogram
X = cat(2, Xs{:});
hist = conHist(X, 'alg', 'eq', 'val', 10);

Cs = zeros(k, k, nS, nS);
for iS = 1 : nS
    for jS = iS + 1 : nS
        Xi = Xs{iS}; segi = segs{iS}; Ii = H2I(segi.st, segi.H);
        Xj = Xs{jS}; segj = segs{jS}; Ij = H2I(segj.st, segj.H);
        
        for ci = 1 : k
            for cj = 1 : k
                Yi = Xi(:, Ii(ci, :) == 1);
                Yj = Xj(:, Ij(cj, :) == 1);
                
                Hi = mapHist(Yi, hist);
                Hj = mapHist(Yj, hist);
                d = conDistHist(Hi, Hj);
        
                Cs(ci, cj, iS, jS) = d;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cs = genConfOver(segs)
% Generate the confusion matrix from multiple segmentations in a fast way.
%
% Input
%   segs  -  segmentation set, 1 x nSub (cell)
%
% Output
%   Cs    -  the confusion matrix, k x k x nSub x nSub

n = length(segs);
k = size(segs{1}.H, 1);
Cs = zeros(k, k, n, n);

for i = 1 : n
    segi = segs{i};
    for j = i + 1 : n
        segj = segs{j};

        C = genConMatrix(segi, segj);
        
        Cs(:, :, i, j) = C;
    end
end
