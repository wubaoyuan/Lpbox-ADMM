function X = gmPosCGrad(K, X0, par)
% Graduate assignment.
%
% Reference
%   S. Gold and A. Rangarajan. "A graduated assignment algorithm
%   for graph matching", IEEE Trans. Pattern Anal. Mach. Intell.,
%   18(4):377–388, 1996.
%
% Math
%   This algorithm is to obtain the optimal X for the following problem
%     max_x   x' * K * x
%     s.t.    A * x <= 1
%
% Remark
%   nn = n1 x n2
%
% Input
%   K        -  affinity matrix, nn x nn (sparse)
%   X0       -  initial assignment, n1 x n2
%   par      -  parameter
%     b0     -  max(n1, n2) | .5
%     bStep  -  {1.075}
%     bMax   -  {200} | 10
%     tolB   -  {1e-3}
%     tolC   -  {1e-3}
%
% Output
%   X        -  permutation matrix, n1 x n2
%
% History
%   create   -  Timothee Cour (timothee.cour@gmail.com), 02-01-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-22-2011

% dimension
[n1, n2] = size(X0);
E12 = ones(n1, n2);

% function parameter
b0 = ps(par, 'b0', max(n1, n2));
bStep = ps(par, 'bStep', 1.075);
bMax = ps(par, 'bMax', 200);
tolB = ps(par, 'tolB', 1e-3);
tolC = ps(par, 'tolC', 1e-3);
maxBIters = 1000;
nthIter = 200;

% affinity
W = tril(K);

% normalize
if 1
    sumW = 2 * sum(sum(W)) - sum(diag(W));
    if sumW == 0
        sumW = eps;
    end
    W = W * length(W) / sumW;
end

indV = find(E12);

% initial solution
X = X0;

% main loop
b = b0;
i = 1;
aIter = 1;
while b < bMax

    oldX{1} = X;
    oldX{2} = X;
    bIter = 1;
    
    % second loop
    while 1
        % X = exp(b * (W * X(:)));
        X(indV) = exp(b * mex_w_times_x_symmetric_tril(X(indV), W));
        X = reshape(X, n1, n2);

        % stop condition
        if any(isnan(X(:)) | isinf(X(:)))
            X = oldX{1};
            break;
        end

        % normalize
        [X, Xslack] = bistocNormalize_slack(X, tolC);

        % stop condition
        if any(isnan(X(:)) | isinf(X(:)))
            X = oldX{1};
            break;
        end

        % prevent oscillations
        % diff = sum(sum(abs(X - oldM)));
        % oldM = X;
        diff1 = sum(abs(X(:) - oldX{1}(:)));
        diff2 = sum(abs(X(:) - oldX{2}(:)));
        diff = min(diff1, diff2);
        oldX{2} = oldX{1};
        oldX{1} = X;

        % stop condition
        if diff < tolB || bIter > maxBIters
            break
        end
        bIter = bIter + 1;
        i = i + 1;
    end
    
    nbIters(aIter) = bIter;
    aIter = aIter + 1;
    b = b * bStep;
end
