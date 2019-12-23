function [X, lambda, timing, constraintViolation] = compute_matching_from_W(W, E12, k, constraintMode, isAffine, isNcut)
% E12 = matching hypothesis (E12 = ones(n1, n2) if full-matching)
%
% Input
%   k = # eigenvectors
%   constraintMode = 'none' | 'col' | 'row' | 'both'
%   isAffine = 0 | 1

W = tril(sparse(W));

constraintViolation = [];
[n1, n2] = size(E12);
n12 = length(W);
% W = tril(W);

if nargin < 6
 isNcut = 0;
end
if nargin < 5
    isAffine = 1;
end

time_eigensolverTotal = cputime;
if strcmp(constraintMode, 'none')
    [X0, lambda, timing] = computeKFirstEigenvectors(W, k, isNcut);

else
    [C, b] = computeConstraintMatching(E12, constraintMode);

    if strcmp(constraintMode, 'both')
        C = C(1 : end - 1, :); % otherwise overconstrained
        b = b(1 : end - 1);
    end

    if isAffine
        [X0, lambda, timing, constraintViolation] = computeEigenvectorsAffineConstraint(W, C, b, k, isNcut);
    else
        [X0, lambda, timing] = computeNcutConstraint_projection(W, C, k, isNcut);
    end
end

% time
time_eigensolverTotal = cputime - time_eigensolverTotal;
timing.eigensolverTotal = time_eigensolverTotal;

X = zeros(n1 * n2, k);
X(E12 > 0, :) = X0;
X = reshape(X, n1, n2, k);
X = X(:, :, 1);
