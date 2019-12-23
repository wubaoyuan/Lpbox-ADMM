function [p, dsts] = meKnl(K)
% Find the centroid in the kernel space.
%
% Input
%   K       -  kernel matrix, n x n
%
% Output
%   p       -  position of centroid
%   dsts    -  distance between each sample and the center, 1 x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 04-22-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

n = size(K, 1);

% sample-center distance
a = sum(K(:)) / n ^ 2;
dsts = zeros(1, n);
for i = 1 : n
    dsts(i) = K(i, i) - 2 * sum(K(i, :)) / n + a;
end

[dst, p] = min(dsts);
