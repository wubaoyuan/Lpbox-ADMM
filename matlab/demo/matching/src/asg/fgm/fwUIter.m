function [X, nIt, objs] = fwUIter(X0, alp, nItMa, isObj)
% Frank-Wolfe iteration for optimizing the symmetric graph matching.
%
% Input
%   X0      -  initial assignment, n1 x n2
%   alp     -  alpha
%   nItMa   -  maximum #iterations
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
global isGXG isHXH;

objs = zeross(1, nItMa);

% main iteration
for nIt = 1 : nItMa
    isGXG = 0;
    isHXH = 0;

    % gradient
    GrVex = fwUGradVex(X0);
    GrCav = fwUGradCav(X0);
    Gr = (1 - alp) * GrVex + alp * GrCav;

    % optimal direction
    YY = gmPosDHun(Gr);
    Y = YY - X0;

    % step size
    [aVex, bVex] = fwUStepVex(X0, Y);
    [aCav, bCav] = fwUStepCav(X0, Y);
    a = (1 - alp) * aVex + alp * aCav;
    b = (1 - alp) * bVex + alp * bCav;
    t = fwStepOpt(a, b);
    
    % update
    X = X0 + t * Y;

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
