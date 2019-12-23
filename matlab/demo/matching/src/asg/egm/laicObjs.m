function [obj, objs, obj1s, Obj2] = laicObjs(X, Hs, Pt, C, bas, lams)
% Compute the objective function of LAIC.
%
% Input
%   X       -  correspondence matrix, Nm x Nt
%   Pt      -  point position, 2 x Nt
%   C       -  feature matching cost matrix, k x n
%   H       -  reconstruction matrix of the template point set, 1 x m (cell), Nm x Nm
%   lams    -  lambda, 1 x k
%
% Output
%   obj     -  objective
%   objs    -  objective for 1st order term, 1 x (k + 1)
%   obj1s   -  k x 1
%   Obj2    -  k x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-13-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 08-28-2012

% dimension
[k, ~] = size(C);
m = length(Hs);
objs = zeros(1, m + 1);

% 1st order term
obj1s = zeros(1, k);
for c = 1 : k
    obj1s(c) = sum(C(c, bas{c}) .* X(c, bas{c})) * lams(c);
end
objs(1) = sum(obj1s);

% 2nd order term
Obj2 = zeros(k, m);
for i = 1 : m
    Obj2(:, i) = sum(abs(Hs{i} * X * Pt'), 2);
    objs(i + 1) = sum(Obj2(:, i));
end

% obj
obj = sum(objs);
