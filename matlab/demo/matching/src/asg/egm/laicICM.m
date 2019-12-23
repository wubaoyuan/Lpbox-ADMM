function X = laicICM(X0, Hs, T, C, bas, lams)
% This function models a linear programming (LP) model for matching and solves it.
%
% Input
%   Hs      -  reconstruction matrix of the template point set, 1 x m (cell), k x k
%   T       -  target points' coordinates, n x d (=2)
%   C       -  feature matching cost matrix, k x n
%   bas     -  candidate matching points, k x 1 (cell), 1 x n
%   lams    -  lambda, 1 x k
%
% Output
%   X       -  correspondence matrix, k x n
%   obj     -  objective
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-13-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 08-28-2012

% dimension
m = length(Hs);
[k, n] = size(C);
d = size(T, 2);

optJs = zeros(1, k);
for c = 1 : k
    bac = bas{c};
    nc = length(bac);
    
    objs = zeros(1, nc);
    for j = 1 : nc
        X = X0;
        X(c, bac) = 0;
        X(c, bac(j)) = 1;
        
        % data term
        obj = sum(sum(C .* X, 2) .* lams(:));
        
        % evaluate the objectve function
        for i = 1 : m
            obj = obj + sum(sum(abs(Hs{i} * X * T), 2));
        end
        
        objs(j) = obj;
    end
    
    [~, optJs(c)] = min(objs);
end

X = zeros(k, n);
for c = 1 : k
    bac = bas{c};
    optJ = optJs(c);
    X(c, bac(optJ)) = 1;
end
