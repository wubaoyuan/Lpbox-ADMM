function Y = mds(D, m, varargin)
% Multi-dimensional Scaling (MDS).
%
% Input
%   D       -  distance matrix, n x n
%   m       -  number of dimensionality after embedding
%   varargin
%     nei   -  neighbour size, {.1}
%
% Output
%   Y       -  new sample matrix after embedding, m x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-26-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

n = size(D, 1);
D2 = -.5 * (D .^ 2 - sum(D .^ 2)' * ones(1, n) / n - ones(n, 1) * sum(D .^ 2) / n + sum(sum(D .^ 2)) / (n ^ 2));

opt.disp = 0;
[vec, val] = eigs(D2, m, 'lm', opt);

Y = vec .* (ones(n, 1) * sqrt(val)')'; 
