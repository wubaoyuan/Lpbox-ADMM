function [X, G] = smpGausses(mes, Vars, ns, varargin)
% Sampling from a set of Gaussian distribution.
%
% Input
%   mes     -  mean, 1 x k (cell), d x 1
%   Vars    -  variance, 1 x k (cell), d x d
%   ns      -  #samples in classes, 1 x k
%
% Output
%   X       -  sample matrix, d x n (= sum(ns))
%   G       -  class indicator matrix, k x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 07-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
k = length(mes);
d = size(mes{1}, 1);

% class size
if length(ns) == 1
    ns = zeros(1, k) + ns;
end

% generate sample
Xs = cell(1, k);
for c = 1 : k
    Xs{c} = smpGauss(mes{c}, Vars{c}, ns(c));
end
X = cat(2, Xs{:});

% label
[~, G] = n2s(ns);
