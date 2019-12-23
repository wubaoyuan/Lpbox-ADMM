function [X, lambda, timing, constraintViolation] = computeEigenvectorsAffineConstraint(W, C, b, nbEigenvectors, isNcut)
%
% Input
%   W       -  sparse upper triangular matrix (n x n)
%              If you have a symmetric matrix W, type W = sparse(tril(W)); before.
%              It will treat W as if it were symmetric, although more efficiently.
%   C       -  constraint matrix (k x n)
%   b       -  constraint matrix (k x 1) (b cannot be all 0)
%   nbEigenvectors  -  #eigenvectors requested
%   isNcut  -  1 or 0 if require normalized cut constrained eigenvectors 
%               instead of constrained eigenvectors
% 
% Output
%   X       -  n x nbEigenvectors (constrained eigenvectors)
%   lambda  -  nbEigenvectors x 1 (eigenvalues)
%   timing  -  computation time
%   constraintViolation  -  error in constraint violation
% 
% solves
%   max sum_k   Xk ^ T W Xk / Xk ^ T Xk
%   subject to  C Xk = b
%               X^T X = I_nbEigenvectors
%   where       Xk = X(:, k)
%   and likewise with normalized cut eigenvectors
%
%   C  -  k x n
%   Ceq = [C(1:k-1,:) - (1/b(k)) * b(1:k-1)*C(k,:)];
%   K = In - Ceq'*inv(Ceq*Ceq')*Ceq;
%   Ceq' * inv(Ceq * Ceq') * Ceq = B' * Ainv * B + L1 * L2'

[C, b] = compute_full_rank_constraint(C, b);

if nargin < 5
    isNcut = 0;
end

W = sparse(W);
C = sparse(C);

if isNcut
    D = mex_computeRowSum(abs(W));
    Dinvsqrt = 1 ./ sqrt(D + eps);
    if issparse(W)
        W = spmtimesd(W, Dinvsqrt, Dinvsqrt);
    else
        W = W .* (Dinvsqrt * Dinvsqrt');
    end
    C = spmtimesd(C, [], Dinvsqrt);
end

% make sure b(end) non zero
[ignore, ind] = max(abs(b));
assert(b(ind) ~= 0);
b([end, ind]) = b([ind, end]);
C([end, ind], :) = C([ind, end], :);

[k, n] = size(C);
B = C(1 : k - 1, :);
u = -b(1 : k - 1) / b(k);
v = C(k, :)';

Ainv = inv(B * B');

T1 = Ainv * [u, B * v + (v' * v) * u];
T2 = T1 * inv(eye(2) + [B * v, u]' * T1);
T3 = -(B' * T2 + v * (u' * T2));

L1 = [B' * (Ainv * u), B' * (Ainv * (B * v)) - v];
L2 = [v + T3(:, 2), T3(:, 1)];

constraint.L1 = full(L1);
constraint.L2 = L2; 
constraint.B = B;
constraint.Ainv = tril(Ainv);

n = size(W, 1);
[options, nbEigenvectors] = getDefaultOptionsEigs(n, nbEigenvectors);
% options.sigma='SA';

options.maxit = 500;
options.tol = 1e-4;
% options.tol=1e-15;

if issparse(W)
    [result, timing] = eigs_compatible_with_eigs_optimized(@mex_projection_affine_symmetric, [], nbEigenvectors, options, tril(W), constraint, []);
%    [result, timing] = eigs_optimized(@mex_projection_affine_symmetric, [], nbEigenvectors, options, tril(W), constraint, []);
else
    error('not implemented');
end
lambda = result.lambda;
X = result.X;

% scaling X to satisfy affine constraint
temp = b' * C * X / (b' * b);
X = X ./ repmat(temp, size(X, 1), 1);
temp = C * X - repmat(b, 1, size(X, 2));
constraintViolation = sqrt(sum(temp .^ 2, 1));

%TODO:do a warning if constraint not satisfied ?
% disp(['constraint violation (before renormalization)= ',num2str(constraintViolation)]);
if isNcut
    X = spdiags(Dinvsqrt, 0, n, n) * X;
end
