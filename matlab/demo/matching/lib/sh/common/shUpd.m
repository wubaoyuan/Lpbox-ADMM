function shUpd(h, X, varargin)
% Show points in 2-D (update).
%
% Input
%   h        -  figure handle
%   X        -  sample matrix, dim x n
%   varargin
%     G      -  class indicator matrix, {[]} | k x n
%
% Output
%   h        -  figure handle
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 12-30-2008
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
G = ps(varargin, 'G', []);

% dimension
X = pca(X, struct('d', 2, 'debg', 'n'));
k = length(h.mks);

% label
if isempty(G)
    G = ones(1, n); 
end
k = size(G, 1); 
l = G2L(G);

% plot
for c = 1 : k
    Y = X(:, l == c);

    set(h.mks{c}, 'XData', Y(1, :), 'YData', Y(2, :));
end
