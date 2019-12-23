function ha = shClu3D(X, G, varargin)
% Show points in clusters in 3-D space.
%
% Remark
%   The dimension of feature, d, has to be 2.
%   To plot points with d > 2, you could call function pca in advance, e.g.,
%      X = pca(X, st('d', 2))
%
% Input
%   X        -  sample matrix, d x n
%   G        -  class indicator matrix, [] | k x n
%   varargin
%     show option
%     parMk  -  marker parameter, {[]}, see function plotmk for more details
%     parAx  -  axis parameter, {[]}, see function setAx for more details
%
% Output
%   ha       -  figure handle
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 12-30-2008
%   modify   -  Feng Zhou (zhfe99@gmail.com), 04-02-2012

% show option
psSh(varargin);

% function option
parMk = ps(varargin, 'parMk', []);
parAx = ps(varargin, 'parAx', []);

% dimension
[d, n] = size(X);
if d ~= 3
    error('unsupported dimension: %d', d);
end

% label matrix (k x n)
if isempty(G)
    G = ones(1, n);
end
k = size(G, 1);
l = G2L(G);

% default marker parameter
if isempty(parMk)
    parMk = st('mkSiz', 5, 'lnWid', 0);
end

% default axis parameter
if isempty(parAx)
    parAx = st('eq', 'y');
end

% main plot
hold on;
ha.mks = cell(1, k);
for c = 1 : k
    ha.mks{c} = plotmk(X(:, l == c), c, parMk);
end

% axis
ha.box = xBox(X, parAx);
setAx(ha.box, parAx);
