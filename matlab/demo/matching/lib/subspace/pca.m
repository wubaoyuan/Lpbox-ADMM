function [Y, V, lam] = pca(X, par)
% Principal Componant Analysis (PCA).
%
% Input
%   X       -  sample matrix, d0 x n
%   par     -  parameter
%     homo  -  flag of using homogeneous in input, 'y' | {'n'}
%     idxD  -  index of dimension, {[]}
%     d     -  parameters for deciding number of dimensions kept after projection, {.999}
%                d < 1:  percentage of energy
%                1 <= d: #dimensions
%     debg  -  debug flag, 'y' | {'n'}
%     fig   -  figure for debug, {1} (only used if debg == 'y')
%
% Output
%   Y       -  principal components, d x n
%   V       -  principal directions, d0 x d
%   lam     -  eigenvalues, d x 1
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-30-2008
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function parameter
isHomo = psY(par, 'homo', 'n');
idxD = ps(par, 'idxD', []);
d = ps(par, 'd', []);
isDebg = psY(par, 'debg', 'n');
fig = ps(par, 'fig', 1);

% input is using homogeneous coordinate
if isHomo
    X(end, :) = [];
end

% used dimension
if ~isempty(idxD)
    X = X(idxD, :);
end

% dimension
[d0, n] = size(X);
if isempty(d) || d0 <= d
    Y = X;
    V = [];
    lam = [];
    return;
end

% centralize
me = sum(X, 2) / n;
X = X - repmat(me, 1, n);

% SVD
[U, S] = svd(X, 0);
[lam0, idx] = sort(sum(S, 2), 'descend');
V0 = U(:, idx);

% energy
if d < 1
    d = thEgy(lam0, d);
end
V = V0(:, 1 : d);
lam = lam0(1 : d);

% projection
Y = V' * X;

% debug
if isDebg
    rows = 1; cols = 1;
    axs = iniAx(fig, rows, cols, [300 * cols, 300 * rows]);
    set(gcf, 'CurrentAxes', axs{1}); hold on;

    plot(1 : length(d0), d0, '.r', 'MarkerSize', 4);
    plot(d, lam0(d), '+b', 'MarkerSize', 7);
    title(sprintf('dim = %d, egy = %.3f', d, sum(lam) / sum(lam0)));
end
