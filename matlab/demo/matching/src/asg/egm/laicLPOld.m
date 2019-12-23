function [X, obj] = laicLPOld(H, T, C, basis, Nsame)
% This function models a linear programming (LP) model for matching and solves it.
%
% Input
%   H       -  reconstruction matrix of the template point set, k x k
%   T       -  target points' coordinates, n x 2
%   C       -  feature matching cost matrix, k x n
%   basis   -  candidate matching points, k x 1 (cell), 1 x ni
%              The ith element records ith model point's candidate matching points in the target point-set
%   Nsame   -  A parameter tuning the "matching to the same target" constraint.
%              For each target point, no more than Nsame model points can be matched to it.
%
% Output
%   X       -  correspondence matrix, k x n
%   obj     -  objective
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-13-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-27-2012

% dimension
[k, n] = size(C);
d = size(T, 2);
lpTool = 'matlab';

% basis
ns = zeros(k, 1);
for i = 1 : k
    ns(i) = length(basis{i});
end
pEds = cumsum(ns);
pHds = [1; pEds(1 : end - 1) + 1];

% formulate matrix for the same point matching constraints
[Vis, Idx] = zeross(k, n);

sumInd = 0;
for i = 1 : k
    Vis(i, basis{i}) = 1;
    Idx(i, basis{i}) = sumInd + 1 : sumInd + ns(i);
    sumInd = sumInd + ns(i);
end

% used target point
nUse = sum(sum(Vis, 1) ~= 0);
idxUse = find(sum(Vis, 1) ~= 0);

% object - basis part
f1 = zeros(1, pEds(end));
for i = 1 : k
    f1(pHds(i) : pEds(i)) = C(i, basis{i});
end
xNum = length(f1);

% object - geometric part
f2 = ones(1, k * d * 2);

% object - combine
f = [f1, f2];
nVar = length(f);

% constraint
A = zeros(k * 3 + nUse, nVar);

% con - sum part (k): sum(X) = 1
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
        T1 = T(basis{j}, 1);
        T2 = T(basis{j}, 2);

        yzInd1 = xNum + i * 4 - 3;
        yzInd2 = yzInd1 + 2;

        A(d1, pHds(j) : pEds(j)) = H(i, j) * T1;
        A(d1, yzInd1 : yzInd1 + 1) = [1 -1];

        A(d2, pHds(j) : pEds(j)) = H(i, j) * T2;
        A(d2, yzInd2 : yzInd2 + 1) = [1 -1];
    end
end

% constraint - mapping multiple points to one point (nUse)
for i = 1 : nUse
     ind = Idx(:, idxUse(i));
     A(k * 3 + i, Idx(ind ~= 0, idxUse(i))) = 1;
end

% b
b = zeros(size(A, 1), 1);
b(1 : k) = 1;
b(k * 3 + 1 : end) = Nsame;

% constraints relations
InEqs = zeros(size(A, 1), 1);
InEqs(k * 3 + 1 : end) = -1;

% solving lp using lp_solver
if strcmp(lpTool, 'lp_maker')
    lp = lp_maker(f, A, b, InEqs, [], [], [], 1, 1);
    solvestat = mxlpsolve('solve', lp);
    obj = mxlpsolve('get_objective', lp);
    compactX = mxlpsolve('get_variables', lp);
    mxlpsolve('delete_lp', lp);

% solving lp using matlab
elseif strcmp(lpTool, 'matlab')
    posEq = 1 : k * 3;
    posIneq = k * 3 + 1 : size(A, 1);

    Aeq = A(posEq, :);
    beq = b(posEq);
    Aineq = A(posIneq, :);
    bineq = b(posIneq);
    lb = zeros(length(f), 1);
    opt.Display = 'off';
    compactX = linprog(f', Aineq, bineq, Aeq, beq, lb, [], [], opt);
    obj = f * compactX;

else
    error('unknown lp toolbox: %s\n', lpTool);
end

% output
X = sparse(k, n);
for i = 1 : k
     X(i, basis{i}) = compactX(pHds(i) : pEds(i));
end
