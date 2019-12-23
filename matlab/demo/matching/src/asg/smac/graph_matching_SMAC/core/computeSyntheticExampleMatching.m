function [W, E12, target, G1, G2] = computeSyntheticExampleMatching(n1, n2, param, noise, density_matches)
% 
% Input
%   n1  -  #nodes in 1st graph
%   n2  -  #nodes in 2nd graph
%   param
%   noise
%   density_matches
% 
% Output
%   W
%   E12
%   target
%   G1
%   G2

% function parameter
if nargin < 3
    param.na = 1; % nb attributes
    param.ratio_fill = .05; % density of edges
end
if nargin < 4
    noise.noise_self = 2; % noise parameter
    noise.noise_edges = 0; % 0.05; %noise structural edge change (sparsity pattern)
    noise.noise_add = 0; % 0.5; 
end
if nargin < 5
    density_matches = 1;
end

% dimension
nMax = max(n1, n2);

% latent graph
[G1, W1] = computeRandomGraph(nMax, param);

% reorder
[G2, W2] = perturbGraph(nMax, G1, W1, noise);
[target, G2, W2] = permuteIndexes(G2, W2, nMax, nMax);

target = target(1 : n1, 1 : n2);
G1 = G1(1 : n1, 1 : n1, :);
G2 = G2(1 : n2, 1 : n2, :);
W1 = W1(1 : n1, 1 : n1);
W2 = W2(1 : n2, 1 : n2);

na = size(G1, 3);
F12 = ones(n1, n2);
E12 = rand(n1, n2) <= density_matches;
E12(target == 1) = 1;
f.fG = ones(na, 1);

options.normalization = 'none';
W = compute_matchingW(G1, W1, G2, W2, F12, E12, f, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [target, G2, W2] = permuteIndexes(G2, W2, n1, n2)
% perm = randperm(n1);
perm = 1 : n1;
G2 = G2(perm, perm, :);
W2 = W2(perm, perm);
target = ind2perm(perm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G1, W1] = computeRandomGraph(n1, param)
%
% Input
%   n1
%   param
%
% Output
%   G1  -  node feature, n1 x n1 x 
%   W1  -  edge feature, n1 x n1 (sparse)

na = param.na;
ratio_fill = param.ratio_fill;

% edge feature
W1 = sprand(n1, n1, ratio_fill) > 0;

% node feature
G1 = computeAttributes(W1, na);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G2, W2] = perturbGraph(n2, G1, W1, noise)
[n1, n1, na] = size(G1);
if n1 ~= n2
    error('TODO');
end

Wnoise = sprand(n2, n2, noise.noise_edges) > 0;
Gnoise_add = computeAttributes(Wnoise, na);
Gnoise_self = computeAttributes(W1, na);
W2 = W1 | Wnoise;
Gnoise = Gnoise_add * noise.noise_add + Gnoise_self * noise.noise_self;
G2 = G1 + Gnoise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G1 = computeAttributes(W1, na)
%
% Input
%   W1  -  edge feature, n1 x n1 (sparse)
%   na  -  
%
% Output
%   G1  -  node feature, n1 x n1 x na

% dimension
n1 = length(W1);
[indi, indj] = find(W1);
ne = length(indi);

for a = 1 : na
    vala = rand(ne, 1);
    %     vala=randn(ne,1);
    G1(:, :, a) = accumarray([indi, indj], vala, [n1, n1]);
end
