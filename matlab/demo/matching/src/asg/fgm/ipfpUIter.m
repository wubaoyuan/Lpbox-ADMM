function [X, nIt] = ipfpUIter(X0, nItMa)
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
global L Ct;
global isHXH;
global IndH2 IndH1T;

isDeb = 1;
if isDeb
    rows = 1; cols = 3;
    Ax = iniAx(9, rows, cols, [250 * rows, 250 * cols]);
    ha = [];
    X00 = X0;
end

% main iteration
objs = zeros(1, nItMa);
for nIt = 1 : nItMa
    isHXH = 0;

    % gradient
    Gr = fwUGradGm(X0);

    % optimal search direction
%    YY = gmPosDHun(Gr, Ct);
    YY = gmPosDHunFast(Gr, Ct);    
    Y = YY - X0;

    % step size
    [a, b] = fwUStepGm(X0, Y);

    % optimal step size
    t = fwStepOpt(a, b);

    % update
    X = X0 + t * Y;
    
    % debug
    if isDeb
        objs(nIt) = multTr(L .* multGXH(IndH1T, X, IndH2) .^ 2);
        ha = deb(ha, Ax, nIt, objs, X, X00);
    end

    % stop condition
    if norm(X(:) - X0(:)) < eps || t < eps    
        break;
    end

    % store
    X0 = X;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ha = deb(ha, Ax, nIt, objs, X, X00)
% Debug
%
% Input

% score
shIt(objs(1 : nIt), ones(1, nIt), 'ax', Ax{1, 1}, 'mkSiz', 7, 'itMa', 0); 
title('objective');

% solution
if nIt == 1
    ha.hX = shM(X, 'ax', Ax{1, 2});
else
    shMUpd(ha.hX, X);
end

% solution
if nIt == 1
    ha.hX00 = shM(X00, 'ax', Ax{1, 3});
end

drawnow;
pause(1);
