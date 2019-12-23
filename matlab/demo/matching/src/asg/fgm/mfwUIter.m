function [X, nIt, objs] = mfwUIter(X0, alp, nItMa, nHst, isObj)
% Modified Frank-Wolfe iteration for optimizing the undirected graph matching.
%
% Input
%   X0      -  initial assignment, n1 x n2
%   alp     -  alpha
%   nItMa   -  maximum #iterations
%   nHst    -  #history node
%   isObj   -  flag of computing objective, 0 | 1
%
% Output
%   X       -  correspondence matrix, n1 x n2
%   nIt     -  #iterations
%   objs    -  objective, 1 x nItMa
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-20-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% global variables
global Ct;
global isGXG isHXH;

objs = zeross(1, nItMa);
Ys = cell(1, nHst);

% main iteration
for nIt = 1 : nItMa
    isGXG = 0;
    isHXH = 0;

    % gradient
    GrVex = fwUGradVex(X0);
    GrCav = fwUGradCav(X0);
    Gr = (1 - alp) * GrVex + alp * GrCav;

    % optimal direction
    Y = gmPosDHun(Gr, Ct);
    V = Y - X0;

    % save to history
    pHst = mod(nIt - 1, nHst) + 1;
    Ys{pHst} = Y / nHst;

    % alternative direction
    if nIt >= nHst
        W = -X0;
        for iHst = 1 : nHst
            W = W + Ys{iHst};
        end

        vV = multTr(Gr .* V) / norm(V, 'fro');
        vW = multTr(Gr .* W) / norm(W, 'fro');
        if vW > vV
            V = W;
            Ys{pHst} = Y / nHst;
        end
    end

    % step size
    [aVex, bVex] = fwUStepVex(X0, V);
    [aCav, bCav] = fwUStepCav(X0, V);
    a = (1 - alp) * aVex + alp * aCav;
    b = (1 - alp) * bVex + alp * bCav;
    t = fwStepOpt(a, b);

    % update
    X = X0 + t * V;

    % debug
    if isObj
        objs(nIt) = pathUObj(X, alp);
    end

    % stop condition
    if norm(X(:) - X0(:)) < eps || t < eps
        break;
    end

    % store
    X0 = X;
end
