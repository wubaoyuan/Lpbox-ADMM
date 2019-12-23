function [Ind, Ind2] = rankX(S, H, m)
% Rank the sample regarding to the distance to the class center.
%
% Input
%   S       -  similarity matrix, n x n
%   H       -  indicator matrix, k x n
%   m       -  number of sample for each class
%
% Output
%   Ind     -  index of samples, n x k
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-25-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

[k, n] = size(H);

Ind = zeros(n, k);
for c = 1 : k
    visC = H(c, :) == 1;
    
    s = sum(S(visC, :), 1);
    [s2, ind] = sort(s, 'descend');
    
    Ind(:, c) = ind';
end

Ind2 = zeros(m, k);
for c = 1 : k
    idx = find(H(c, :));
    nc = length(idx);
    
    vec = Ind(1 : nc, c);
    blocks = vec2block(vec, m - 1);
    
    for i = 1 : m - 1
        Ind2(i, c) = blocks{i}(1);
    end
    Ind2(m, c) = blocks{end}(end);
end
