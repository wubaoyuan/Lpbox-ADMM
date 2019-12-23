function [KP, KQ, KS, KA] = conKnlGph(XPs, XQs, gphs, parKnl)
% Compute affinity matrix for graph matching.
%
% Remark
%   nnT = n1T x n2T
%   If visPts{1} == [], then visPts{1} is default to ones(ki, 1)
%   If visPts{1} ~= [], then ki = length(find(visPts{i}))
%
% Input
%   XPs       -  node feature, 1 x 2 (cell), dP x ni
%   XQs       -  edge feature, 1 x 2 (cell), dQ x 2mi
%   gphs      -  graphs, 1 x 2 (cell)
%   parKnl    -  parameter
%     visPts  -  node existence status, {[]} | 1 x 2 (cell), niT x 1
%     alg     -  method of computing affinity, {'toy'} | 'cmum' | 'pas'
%                'toy':  toy data
%                'cmum': CMU motion data
%                'pas':  Pascal data
%
% Output
%   KP        -  node-node affinity, n1T x n2T
%   KQ        -  edge-edge affinity, m1 x m2
%   KS        -  affinity (based on symmetric edges), nnT x nnT (sparse)
%   KA        -  affinity (based on asymmetric edges), nnT x nnT (sparse)
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 08-09-2011
%   modify    -  Feng Zhou (zhfe99@gmail.com), 01-29-2012

% function parameter
visPts = ps(parKnl, 'visPts', []);
alg = ps(parKnl, 'alg', 'toy');

% dimension
n1 = size(XPs{1}, 2);
n2 = size(XPs{2}, 2);
m1 = size(XQs{1}, 2);
m2 = size(XQs{2}, 2);
prIn('conKnlGph', 'alg %s, n1 %d, n2 %d, m1 %d, m2 %d', alg, n1, n2, m1, m2);

% node & edge affinity
[KP, KQ] = knlPQ(XPs, XQs, Egs, alg);

% normalize
% KQ = knlEgNor(KQ, parKnl);

% fill-in the missing data
[KP, Eg2s] = knlFill(KP, Egs, visPts);

% convert edge affinity to node-pair affinity
if nargout == 3
    K = knlPQ2K(KP, KQ, Eg2s);
end

% for my method
%KQ = KQ(1 : round(m1 / 2), 1 : round(m2 / 2));

% normalize
% K = knlEgNor2(K, [n1 n2]);

prOut;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KP, KQ] = knlPQ(XPs, XQs, Egs, alg)
% Compute node and edge affinity.
%
% Input
%   XPs  -  node feature, 1 x 2 (cell), dP x ni
%   XQs  -  edge feature, 1 x 2 (cell), dQ x mi
%   Egs  -  graph edges, 1 x 2 (cell), 2 x mi
%   alg  -  method of computing the kernel, 'toy' | 'cmum' | 'pas'
%           'toy':  toy data
%           'cmum': cmu motion sequence
%           'pas':  Pascal image data
%
% Output
%   KP   -  node affinity, n1 x n2
%   KQ   -  edge affinity, m1 x m2

% dimension
n1 = size(XPs{1}, 2);
n2 = size(XPs{2}, 2);
m1 = size(XQs{1}, 2);
m2 = size(XQs{2}, 2);

% for toy data
if strcmp(alg, 'toy')
    KP = zeros(n1, n2);

    DQ = conDst(XQs{1}(1 : 2, :), XQs{2}(1 : 2, :));
    KQ = exp(-DQ / .15);

% for shp1 data
elseif strcmp(alg, 'shp1')
    KP = zeros(n1, n2);
    
    X1 = XQs{1}(3 : end, :);
    X2 = XQs{2}(3 : end, :);
    X1 = X1 ./ repmat(sum(X1 .^ 2, 1), 2, 1);
    X2 = X2 ./ repmat(sum(X2 .^ 2, 1), 2, 1);    

    % uniform
    KQ = X1' * X2;
    
% for shp2 data
elseif strcmp(alg, 'shp2')
    KP = zeros(n1, n2);
    
    % distance
    DQ = conDst(XQs{1}, XQs{2});
    KQ = exp(-DQ / 2500);
    
% for shp3 data
elseif strcmp(alg, 'shp3')
    KP = zeros(n1, n2);
    
    % angle
    Ang1 = repmat(XQs{1}(2, :)', 1, m2);
    Ang2 = repmat(XQs{2}(2, :), m1, 1);
    Ang = abs(Ang1 - Ang2);
    KQ = exp(-Ang);

% for cmu motion data
elseif strcmp(alg, 'cmum')
    KP = zeros(n1, n2);

    DQ = conDst(XQs{1}, XQs{2});
    KQ = exp(-DQ / 2500);
    
% for cmum
elseif strcmp(alg, 'cmum2')
    KP = zeros(n1, n2);

    % distance
    Dst1 = repmat(XQs{1}(1, :)', 1, m2);
    Dst2 = repmat(XQs{2}(1, :), m1, 1);
    Dst = abs(Dst1 - Dst2) ./ (min(Dst1, Dst2) + eps);

    % angle
    Ang1 = repmat(XQs{1}(2, :)', 1, m2);
    Ang2 = repmat(XQs{2}(2, :), m1, 1);
    Ang = abs(Ang1 - Ang2);
    
    % combine distance and angle
    KQ = exp(-(Dst + Ang) / 2);

% for Pascal data    
elseif strcmp(alg, 'pas')
    KP = zeros(n1, n2);
    
    % contour normal
    XP1 = XPs{1};
    XP2 = XPs{2};
    Eg1 = Egs{1};
    Eg2 = Egs{2};
    I11 = repmat(Eg1(1, :)', 1, m2);
    I12 = repmat(Eg1(2, :)', 1, m2);
    I21 = repmat(Eg2(1, :), m1, 1);
    I22 = repmat(Eg2(2, :), m1, 1);
    Nor1 = abs(XP1(I11) - XP2(I21));
    Nor2 = abs(XP1(I12) - XP2(I22));
    
    % distance
    Dst1 = repmat(XQs{1}(1, :)', 1, m2);
    Dst2 = repmat(XQs{2}(1, :), m1, 1);
    Dst = abs(Dst1 - Dst2) ./ (min(Dst1, Dst2) + eps);

    % angle
    Ang1 = repmat(XQs{1}(2, :)', 1, m2);
    Ang2 = repmat(XQs{2}(2, :), m1, 1);
    Ang = abs(Ang1 - Ang2);
    
    % combine distance and angle
    KQ = exp(-(Nor1 + Nor2 + Dst + Ang) / 4);
    
% for Pascal data    
elseif strcmp(alg, 'pas2')
    KP = zeros(n1, n2);
    
    % distance
    Dst1 = repmat(XQs{1}(1, :)', 1, m2);
    Dst2 = repmat(XQs{2}(1, :), m1, 1);
    Dst = abs(Dst1 - Dst2) ./ (min(Dst1, Dst2) + eps);

    % angle
    Ang1 = repmat(XQs{1}(2, :)', 1, m2);
    Ang2 = repmat(XQs{2}(2, :), m1, 1);
    Ang = abs(Ang1 - Ang2);
    
    % combine distance and angle
    KQ = exp(-(Dst + Ang) / 2);
    
% for Pascal data    
elseif strcmp(alg, 'pas3')
    DP = conDst(XPs{1}, XPs{2});
    KP = exp(-real(sqrt(DP)));
    
    % distance
    Dst1 = repmat(XQs{1}(1, :)', 1, m2);
    Dst2 = repmat(XQs{2}(1, :), m1, 1);
    Dst = abs(Dst1 - Dst2) ./ (min(Dst1, Dst2) + eps);

    % angle
    Ang1 = repmat(XQs{1}(2, :)', 1, m2);
    Ang2 = repmat(XQs{2}(2, :), m1, 1);
    Ang = abs(Ang1 - Ang2);
    
    % combine distance and angle
    KQ = exp(-(Dst + Ang) / 2);
    
else
    error('unknown algorithm: %s', alg);
end
