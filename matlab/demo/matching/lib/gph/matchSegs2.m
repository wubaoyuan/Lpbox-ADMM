function [segs, idxs] = matchSegs2(seg0s, Ss, varargin)
% Find the optimum matching among segmentations.
% Notice:
%   The class number of each segmentation must be the same.
%
% Input
%   Xs      -  set of sample matrix, 1 x nS (cell), dim x nF
%   seg0s   -  set of segmentation, 1 x nS (cell)
%   varargin
%     conf  -  method to create confusion measure between segments, 'dtak' | 'space' | 'over'
%              'dtak':  Dynamic Time Alignment Kernel
%              'space': Spatial Distance
%              'peak':  Peak frame
%              'hist':  Histogram
%              'over':  Overlapping of different segmentations on the same sequence. In this case, Xs is useless.
%     val   -  value for confusion matrix
%     alg   -  algorithm for finding the matching, see function assignM for more details.
%
% Output
%   segs    -  segmentation after matching, 1 x nS (cell)
%   C0s     -  confusion matrix before adjusting, k x k x nS x nS
%   Cs      -  confusion matrix after adjusting, k x k x nS x nS
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-18-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

prompt('matching segments...\n');

% confusion matrix
C0s = genConf(Ss, seg0s);

% multidimensional assignment
[Ps, Cs] = assignM(C0s);

% adjust segmentation
nS = length(seg0s);
H0s = cellFld(seg0s, 'H');
ks = cellDim(H0s, 1);
ms = cellDim(H0s, 2);
[mak, ind] = max(ks);

segs = cell(1, nS);
idxs = cell(1, nS);
for i = 1 : nS
    if i == ind
        segs{i} = seg0s{i};
        continue;
    end

    P = Ps(:, :, ind, i);

    H1 = zeros(mak, ms(i));
    H1(1 : ks(i), :) = H0s{i};

    idx = zeros(1, mak);
    for c = 1 : mak
        p = P(c, :);
        idx(c) = find(p);
    end
    H = H1(idx, :);
    idxs{i} = idx;

    segs{i} = seg0s{i};
    segs{i}.H = H;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cs = genConf(Ss, segs)
% Compute the class confusion matrix for pair of segmentation based on the temporal similarity.
%
% Input
%   Xs      -  set of sample matrix, 1 x nS (cell), dim x nF
%   segs    -  set of segmentation, 1 x nS (cell)
%
% Output
%   Cs      -  confusion matrix, mak x mak x nS x nS
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-29-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

nS = length(segs);
Hs = cellFld(segs, 'H');
ks = cellDim(Hs, 1);
mak = max(ks);

Cs = zeros(mak, mak, nS, nS);
for iS = 1 : nS
    Hi = Hs{iS}; ki = ks(iS); segi = segs{iS};
    Cs(1 : ki, 1 : ki, iS, iS) = Hi * Hi';

    for jS = iS + 1 : nS
        Hj = Hs{jS}; kj = ks(jS); segj = segs{jS};

        % segment similarity
        S = Ss{iS, jS};
        sS = conSegSim(S, segi, segj);

        % sum of pairwise similarity
        C = Hi * sS * Hj';
        Cs(1 : ki, 1 : kj, iS, jS) = C;
        
        Cs(:, :, jS, iS) = Cs(:, :, iS, jS)';
    end
end
