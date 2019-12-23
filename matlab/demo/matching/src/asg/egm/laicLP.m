function [X, obj, res1, res2, A, b] = laicLP(H, T, C, bas, lam)
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
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-31-2012

% dimension
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
    f1(pHds(i) : pEds(i)) = C(i, bas{i});
end
xNum = length(f1);

% object - geometric part
f2 = zeros(1, k * d * 2) + lam;

% object - combine
f = [f1, f2];
nVar = length(f);

% constraint
A = zeros(k + k * d, nVar);

% constraint - sum part: sum(X) = 1
for i = 1 : k
    A(i, pHds(i) : pEds(i)) = 1;
end

% constraint - absolute linearization part (2k)
for i = 1 : k
    % row number
    d1 = k + i * 2 - 1;
    d2 = d1 + 1;

    % for each item in one row of H
    for j = 1 : k
        T1 = T(bas{j}, 1);
        T2 = T(bas{j}, 2);

        yzInd1 = xNum + i * 4 - 3;
        yzInd2 = yzInd1 + 2;

        A(d1, pHds(j) : pEds(j)) = H(i, j) * T1;
        A(d1, yzInd1 : yzInd1 + 1) = [1 -1];

        A(d2, pHds(j) : pEds(j)) = H(i, j) * T2;
        A(d2, yzInd2 : yzInd2 + 1) = [1 -1];
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
    opt.Display = 'off';
    res = linprog(f', [], [], A, b, lb, [], [], opt);
    obj = f * res;
    
% solving lp using CVX
elseif strcmp(lpTool, 'cvx')
    lb = zeros(length(f), 1);    
    cvx_begin
        variable x(nVar, 1)
        minimize(f * x)
        subject to
            A * x == b
            x >= lb
    cvx_end    
    res = x;
    obj = f * res;
else
    error('unknown lp toolbox: %s\n', lpTool);
end

res1 = res(1 : pEds(end));
res2 = res(pEds(end) + 1 : end);

% output
X = zeros(k, n);
for i = 1 : k
     X(i, bas{i}) = res(pHds(i) : pEds(i));
end
XPos = reshape(res2(1 : k * d), [k d]);
XNeg = reshape(res2(k * d + 1 : end), [k d]);

equal('XPos + XNeg', XPos + XNeg, H * X * T);

% test
[res1a, res2a] = laicLPX2Res(X, H, T, bas);
equal('res1', res1, res1a);
equal('res2', res2, res2a);
equal('b', A * [res1a; res2a], b);