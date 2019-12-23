function [X, nIt] = ipfpDIter(X0, nItMa)
% IPFP iteration for optimizing the asymmetric graph matching.
%
% Input
%   X0      -  initial assignment, n1 x n2
%   nItMa   -  maximum #iterations
%
% Output
%   X       -  correspondence matrix, n1 x n2
%   nIt     -  #iterations
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-20-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% global variables
global isGXG isHXH;

% main iteration
for nIt = 1 : nItMa
    isGXG = 0;
    isHXH = 0;

    % gradient
    Gr = fwDGradGm(X0);

    % optimal search direction
    YY = gmPosDHun(Gr);
    Y = YY - X0;

    % step size
    [a, b] = fwDStepGm(X0, Y);

    % optimal step size
    t = fwStepOpt(a, b);

    % update
    X = X0 + t * Y;
    
    % stop condition
    if norm(X(:) - X0(:)) < eps || t < eps    
        break;
    end
    
    % store
    X0 = X;
end
