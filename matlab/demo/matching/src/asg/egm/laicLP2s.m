function [X, obj] = laicLP2s(H, T, C, bas, lam, weis)
% This function models a linear programming (LP) model for matching and solves it.
%
% Input
%   H       -  reconstruction matrix of the template point set, k x k
%   T       -  target points' coordinates, n x d (=2)
%   C       -  feature matching cost matrix, k x n
%   bas     -  candidate matching points, k x 1 (cell), 1 x ni
%   lam     -  lambda, 1 x 1
%   weis    -  weight, 1 x k
%
% Output
%   X       -  correspondence matrix, k x n
%   obj     -  objective
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-13-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 08-02-2012

% dimension
m = 1;
[k, n] = size(C);
d = size(T, 2);
lpTool = 'matlab';

% basis
ns = zeros(k, 1);
for i = 1 : k
    ns(i) = length(bas{i});
end
pEds = cumsum(ns);
pHds = [1; pEds(1 : end - 1) + 1];

% object - basis part
f1 = zeros(1, pEds(end));
for i = 1 : k
    f1(pHds(i) : pEds(i)) = C(i, bas{i}) * lam;
end
xNum = length(f1);

% object - geometric part
f2 = zeros(1, k * d * 2 * m) + 1;

% object - combine
f = [f1, f2];
nVar = length(f);

% constraint
A = zeros(k + k * d * m, nVar);

% constraint - sum part: sum(X) = 1
for i = 1 : k
    A(i, pHds(i) : pEds(i)) = 1;
end

% constraint - absolute linearization part (2km)
[IdxPs, IdxNs] = cellss(1, m);
for c = 1 : m
    [IdxP, IdxN] = zeross(k, 2);
    
    for i = 1 : k
        % row number
        d1 = k + d * k * (c - 1) + i * 2 - 1;
        d2 = d1 + 1;

        % for each item in one row of H
        for j = 1 : k
            T1 = T(bas{j}, 1);
            T2 = T(bas{j}, 2);
            A(d1, pHds(j) : pEds(j)) = Hs{c}(i, j) * T1;
            A(d2, pHds(j) : pEds(j)) = Hs{c}(i, j) * T2;
        end
        yzInd1 = xNum + d * 2 * k * (c - 1) + i * 4 - 3;
        yzInd2 = yzInd1 + 2;

        A(d1, yzInd1 : yzInd1 + 1) = [1 -1];
        A(d2, yzInd2 : yzInd2 + 1) = [1 -1];
        
        IdxP(i, :) = [yzInd1, yzInd2];
        IdxN(i, :) = [yzInd1, yzInd2] + 1;
    end

    IdxPs{c} = IdxP;
    IdxNs{c} = IdxN;
    for i = 1 : k
        f(IdxP(i, :)) = f(IdxP(i, :)) * weis(i);
        f(IdxN(i, :)) = f(IdxN(i, :)) * weis(i);
    end
end

% b
b = zeros(size(A, 1), 1);
b(1 : k) = 1;

% solving lp using lp_solver
if strcmp(lpTool, 'lp_maker')
    lp = lp_maker(f, A, b, InEqs, [], [], [], 1, 1);
    solvestat = mxlpsolve('solve', lp);
    obj = mxlpsolve('get_objective', lp);
    res = mxlpsolve('get_variables', lp);
    mxlpsolve('delete_lp', lp);

% solving lp using matlab
elseif strcmp(lpTool, 'matlab')
    lb = zeros(length(f), 1);
    options = optimset('LargeScale', 'off', 'Simplex', 'on', 'Display', 'off');
    res = linprog(f', [], [], A, b, lb, [], [], options);
    obj = f * res;

else
    error('unknown lp toolbox: %s\n', lpTool);
end

% output
X = zeros(k, n);
for i = 1 : k
     X(i, bas{i}) = res(pHds(i) : pEds(i));
end
[XPs, XNs] = cellss(1, m);
for i = 1 : m
    [XP, XN] = zeross(k, 2);
    IdxP = IdxPs{i};
    IdxN = IdxNs{i};
    for c = 1 : k
        XP(c, :) = res(IdxP(c, :));
        XN(c, :) = res(IdxN(c, :));
    end
    XPs{i} = XP;
    XNs{i} = XN;
end

% testing
% obj2 = sum(sum(C .* X * lam));
% for i = 1 : m
%     obj2 = obj2 + sum(sum(XPs{i} + XNs{i}, 2) .* weis(:));
% end
% equal('obj', obj, obj2);