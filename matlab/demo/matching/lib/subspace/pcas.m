function Ys = pcas(Xs, par)
% Principal Componant Analysis (PCA) for multiple sets of points.
%
% Input
%   Xs      -  set of sample matrices, 1 x m (cell), dim0 x ni
%   par     -  parameter, see function pca for the details
%     cat   -  flag of concatenating X, {'y'} | 'n'
%
% Output
%   Ys      -  set of principal components, 1 x m (cell), dim x ni
%   V       -  principal directions, dim0 x dim
%   lam     -  eigenvalues, dim x 1
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 09-12-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function parameter
isCat = psY(par, 'cat', 'y');

% dimension
m = length(Xs);
Ys = cell(1, m);

if isCat
    ns = cellDim(Xs, 2);
    s = n2s(ns);

    % concatenate
    X = cat(2, Xs{:});

    % pca
    Y = pca(X, par);

    % split
    for i = 1 : m
        Ys{i} = Y(:, s(i) : s(i + 1) - 1);
    end
else
    for i = 1 : m
        Ys{i} = pca(Xs{i}, par);
    end
end
