function Xorth = computeXorthonormal(X, E12)

[n1, n2] = size(X);
[U, S, V] = svd(X);
Xorth = U * eye(n1, n2) * V';
Xorth(E12 == 0) = 0;
