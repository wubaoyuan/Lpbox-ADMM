function X = genCentre(me, var, n, minDist)
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

if ~exist('minDist', 'var'), minDist = 0; end

maxIterN = 30;
dim = size(me, 1);
X = zeros(dim, n);

for i = 1 : n

    joinCoMin = i;
    for iter = 1 : maxIterN
        xi = genGauss(me, var, 1);

        joinCo = 0;
        for j = 1 : i - 1
            xj = X(:, j);
            dist = sqrt((xi - xj)' * (xi - xj));

            if dist < minDist, joinCo = joinCo + 1; end
        end

        if joinCo < joinCoMin
            X(:, i) = xi;
            joinCoMin = joinCo;
        end
    end
end
