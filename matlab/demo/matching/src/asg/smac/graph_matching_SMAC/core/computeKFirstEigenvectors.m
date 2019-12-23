function [X, lambda, timing, D] = computeKFirstEigenvectors(W, nbEigenvectors, isncut, D)
% Computes first k (normlalized-cuts) eigenvectors of W.
%
% Input
%   W  -  nxn affinity matrix
%   nbEigenvectors: # eigenvectors requested (k)
%   isncut: 0 for eigenvectors, 1 for normlalized-cuts eigenvectors
%   D: optional parameter, can be set to diagonal of degree matrix of W; it is computed when not specified
%   X: nxk eigenvectors
%   lambda: kx1 eigenvalues
%   timing: timing information
%
% Output

if nargin < 3
    isncut = 1;
end

% dimension
n = size(W, 1);
[options, nbEigenvectors] = getDefaultOptionsEigs(n, nbEigenvectors);

if isncut
    if nargin >= 4
        D = D(:);
        Dinvsqrt = 1 ./ sqrt(D + eps);
        W = normalizeW_D(W, Dinvsqrt);
    else
        [W, Dinvsqrt, D] = normalizeW_D(W, [], 0);
    end
end

% core algorithm
if issparse(W) && mex_istril(W)    
    [result, timing] = eigs_optimized(@mex_w_times_x_symmetric_tril, [], nbEigenvectors, options, W);
else
    [result, timing] = eigs_optimized(W, [], nbEigenvectors, options);
end
X = result.X;

if isncut
    X = spdiag(Dinvsqrt) * X;
end

lambda = result.lambda;
mex_normalizeColumns(X);
