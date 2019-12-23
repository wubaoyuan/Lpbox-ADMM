function K = cenK(K)
% Centralize the kernel matrix.
%
% Input
%   K       -  original kernel matrix, n x n
%
% Output
%   K       -  new kernel matrix, n x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-26-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 04-20-2013

% dimension
n = size(K, 1);

% matrix form:
%   P = eye(n) - ones(n, n) / n;
%   K = P * K * P;
vals = sum(K, 1);
for i = 1 : n
    K(:, i) = K(:, i) - vals(i) / n;
end
for i = 1 : n
    K(i, :) = K(i, :) - vals(i) / n;
end
K = K + sum(vals) / n ^ 2;
