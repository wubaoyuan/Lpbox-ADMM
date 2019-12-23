function X = norX(X0, nei)
% Normalize each dimension spearately with respect to the variance of distance.
%
% Input
%   X0      -  original sample matrix, dim x n
%   nei     -  width of neighborhood, dim x 1, default to ones(dim, 1)
%
% Output
%   X       -  new sample matrix, dim x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

[dim, n] = size(X0);

% variance in each dimension
vs = ones(dim, 1);
for di = 1 : dim
    % Estimate the variance of distance
    m = min(n - 1, floor(n * nei(di)));
    
    if m > 0
        D = conDist(X0(di, :), X0(di, :), 'dst', 'e');
        D1 = sort(D);
        D2 = D1(2 : (1 + m), :);
        
        vs(di) = sum(sum(D2, 1)) / (n * m);
    end
end
X = X0 ./ repmat(vs, 1, n);
