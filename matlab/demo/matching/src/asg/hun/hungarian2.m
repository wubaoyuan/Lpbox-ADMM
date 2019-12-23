function [X, score] = hungarian2(A);
assert(isa(A, 'double'));
assert(~issparse(A));

[X, score] = assignmentoptimal(A);
[n1, n2] = size(A);
ind1 = find(X);
ind2 = X(ind1);
X = full(sparse(ind1, ind2, 1, n1, n2));
