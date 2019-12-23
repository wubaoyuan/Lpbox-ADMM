function [A, Dire] = rca(X0, H)
% Relevant Component Analysis (RCA).
%
% Input
%   X0      -  sample matrix, dim x n
%   H       -  indicator matrix, k x n
%
% Output
%   A       -  the RCA suggested Mahalanobis matrix, dim x dim
%   Dire    -  direction, dim x dim
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-26-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

[dim, n] = size(X0);
k = size(H, 1);

% subtract the mean
X = normX(X0);

% within-class scatter
sw = zeros(dim);
for c = 1 : k
    Y = X(:, logical(H(c, :)));
    
    Y = normX(Y);

    sw = sw + Y * Y';
end
sw = sw / n;

Dire = sw ^ (-.5);

A = Dire * Dire';
