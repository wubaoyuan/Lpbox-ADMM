function [obj, objs, obj1s, Obj2] = laicObjClus(X, Pt2, C, Hs, wClus, lam, weis)
% Compute the objective function of LAIC.
%
% Input
%   X       -  correspondence matrix, Nm x Nt
%   Pt2     -  point position, 2 x Nt
%   C       -  feature matching cost matrix, k x n
%   H       -  reconstruction matrix of the template point set, 1 x m (cell), Nm x Nm
%   lam     -  lambda, 1 x 1
%   weis    -  weight, 1 x k
%
% Output
%   obj     -  objective
%   objs    -  objective for 1st order term, 1 x (k + 1)
%   obj1s   -  k x 1
%   Obj2    -  k x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-13-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 08-02-2012

% dimension
[k, ~] = size(C);
m = length(Hs);
objs = zeros(1, m + 1);

% 1st order term
obj1s = sum(C .* X, 2) * lam;
objs(1) = sum(obj1s);

% 2nd order term
Obj2 = zeros(k, m);
for i = 1 : m
    if isempty(Hs{i})
        continue;
    end

    Obj2(:, i) = sum(abs(Hs{i} * X * Pt2'), 2) .* weis(:) * wClus(i);
    objs(i + 1) = sum(Obj2(:, i));
end

% obj
obj = sum(objs);
