function [result, timing] = eigs_compatible_with_eigs_optimized(A, B, k, options, Afun1, Afun2, Afun3);

sigma = 'LA';
opts.issym = 1;
opts.disp = 0; 

[opts,k] = getDefaultOptionsEigs(options.n,k);

if isnumeric(A)
    [X, lambda] = eigs(A, B, k, sigma, opts);
else
%    save('bad.mat', 'A', 'options', 'B', 'k', 'sigma', 'opts', 'Afun1', 'Afun2', 'Afun3');
    [X, lambda] = eigs(A, options.n, B, k, sigma, opts, Afun1, Afun2, Afun3);
end

lambda = diag(lambda);

if strcmp(sigma, 'LA')
    [val, indexesOrdered] = sort(lambda,'descend');
elseif strcmp(sigma, 'SA')
    [val, indexesOrdered] = sort(lambda,'ascend');
else
    assert(0);
end
lambda = lambda(indexesOrdered);
X = X(:, indexesOrdered);
result.X = X;
result.lambda = lambda;

timing = [];
% n                  :  length of matrix represented by A if A is a function
% tol                :  default eps
% p                  :  default min(max(2*k,20),n)
% maxit              :  default max(300,ceil(2*n/max(p,1)))
% v0                 :  default zeros(n,1)
% cholB              :  [{0} | 1]
% permB              :  default 1:n
%
% output :
% result.X
% result.lambda
% result.flag
% result.nbA_times_X
% result.nbIterations
% result.nbEigenvaluesConverged
%
% timing.preprocessing
% timing.dsaupd
% timing.A_times_X
% timing.postprocessing
% timing.total
