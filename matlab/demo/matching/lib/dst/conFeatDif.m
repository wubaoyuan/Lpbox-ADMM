function [W, M] = conFeatDif(X, H, varargin)
% Construct the feature difference matrix.
%
% Input
%   X       -  sample matrix, dim x n
%   H       -  indicator matrix, k x n
%   varargin
%     type  -  type, {'a'} | 's' | 'd'
%              'a': all samples
%              's': only similar samples
%              'd': only dis-similar samples
%
% Output
%   W
%   M
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-18-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
type = ps(varargin, 'type', 'a');

[dim, n] = size(X);
W = zeros(dim, n * n);

% I0(i, j) == 1 iff X(:, i) and X(:, j) have the same label
I0 = H' * H;
if strcmp(type, 'a')
    I = ones(size(I0)) == 1;
elseif strcmp(type, 's')
    I = I0 == 1;
elseif strcmp(type, 'd')
    I = I0 == 0;
else
    error('unknown type');
end
    
m = 0; M = zeros(dim, dim);
for i = 1 : n
    for j = 1 : n

        if I(i, j)
            m = m + 1;
            d = X(:, i) - X(: ,j);
            W(:, m) = d;
            M = M + d * d';
        end
    end
end
W(:, m + 1 : end) = [];

