function [Xs, mi, ma] = scaX(X0s, sca)
% Re-scale each dimension of X to [0, sca]
%
% Input
%   X0      -  original sample matrix, dim x n
%   sca     -  numeric range of each dimensions, 1 x 2
%
% Output
%   X       -  new sample matrix, dim x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-23-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
m = length(X0s);
dim = size(X0s{1}, 1);
ns = cellDim(X0s, 2);

% boundary
mi = inf(dim, 1);
ma = -inf(dim, 1);
for i = 1 : m
    X0 = X0s{i};
    mi = min(min(X0, [], 2), mi);
    ma = max(max(X0, [], 2), ma);
end

Xs = cell(1, m);
for i = 1 : m
    Xs{i} = sca * (X0s{i} - repmat(mi, 1, ns(i))) ./ (repmat(ma - mi, 1, ns(i)) + eps);
end
