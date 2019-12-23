function [ind, w] = qda(X, H, lambda)
% Quadratic discriminant analysis (QDA).
%
% Input
%   X       -  sample matrix, dim x n
%   H       -  indicator matrix, k x n
%   lambda  -  lambda, {1}
%
% Output
%   ind     -  feature index (sorted descendly), dim x 1
%   w       -  weight (sorted descendly), dim x 1
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

if ~exist('lambda', 'var')
    lambda = 1;
end

[dim, n] = size(X);
k = size(H, 1);

% move center to zero
mean = sum(X, 2) / n;
X = X - repmat(mean, 1, n);

% total within scatter
st = X * X';

% within-class scatter
sw = zeros(dim);
for c = 1 : k
    Y = X(:, logical(H(c, :)));
    
    m = size(Y, 2);
    mean = sum(Y, 2) / m;
    Y = Y - repmat(mean, 1, m);

    sw = sw + Y * Y';
end

% between-class scatter
sb = st - sw;

% convex programming
H2 = sw - lambda * sb;
Aeq = ones(1, dim); beq = 1;
lb = zeros(dim, 1); ub = ones(dim, 1);
b = quadprog(H2, [], [], [], Aeq, beq, lb, ub);

[w, ind] = sort(b, 'descend');
