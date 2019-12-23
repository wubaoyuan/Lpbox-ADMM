function vars = genVar(dim, n, varargin)
% Generate several variances with the specific parameters.
%
% Input
%   dim     -  dimensionality
%   n       -  #samples
%   varargin
%     mi    -  minimum radius in each dimension, {0}
%     ma    -  maximum radius in each dimension, {1}
%     diag  -  diagonal flag, {'n'} | 'y'
%
% Output
%   vars    -  variances, dim x dim x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-28-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
mi = ps(varargin, 'mi', 0);
ma = ps(varargin, 'ma', 1);
isDiag = psY(varargin, 'diag', 'n');

vars = zeros(dim, dim, n);

for i = 1 : n
    if isDiag
        d = rand(dim, 1);
        x = diag(mi + (ma - mi) * d);
    end
%     x = x + diag(ones(dim, 1) + mi);
%     x = x * ave;

%     vars(:, :, i) = x' * x;
    vars(:, :, i) = x;
end
