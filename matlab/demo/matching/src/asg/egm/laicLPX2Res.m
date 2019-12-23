function [res1, res2, XPos, XNeg] = laicLPX2Res(X, H, T, bas)
% This function models a linear programming (LP) model for matching and solves it.
%
% Input
%   H       -  reconstruction matrix of the template point set, k x k
%   T       -  target points' coordinates, n x d (=2)
%   C       -  feature matching cost matrix, k x n
%   bas     -  candidate matching points, k x 1 (cell), 1 x ni
%   lam     -  lambda
%
% Output
%   X       -  correspondence matrix, k x n
%   obj     -  objective
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-13-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-27-2012

% dimension
[k, n] = size(X);
d = size(T, 2);

Dif = H * X * T;

XPos = (abs(Dif) + Dif) / 2;
XNeg = (abs(Dif) - Dif) / 2;

res1 = [];
for i = 1 : k
    res1 = [res1; X(i, bas{i})'];
end

res2 = [XPos(:); XNeg(:)];
