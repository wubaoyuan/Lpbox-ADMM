function [X, nIt, objs] = mfwDIter(X0, alp, nItMa, nHst, isObj)
% Modified Frank-Wolfe iteration for optimizing the directed graph matching.
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
global isGXG isHXH;

objs = zeross(1, nItMa);
Ys = cell(1, nHst);

X0 = sparse(X0);

% main iteration
for nIt = 1 : nItMa
    
%    fprintf('sparesness of X %.2f\n', nnz(X0)/prod(size(X0)));

    isGXG = 0;
    isHXH = 0;

    % gradient
    GrGm = fwDGradGm(X0);
    GrCon = fwDGradCon(X0);
    Gr = GrGm + (alp - .5) * GrCon;

    % optimal direction
    Y = gmPosDHun(Gr);
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
    
%    fprintf('sparesness of V %.2f\n', nnz(V)/prod(size(V)));

    % step size
    [aGm, bGm] = fwDStepGm(X0, V);
    [aCon, bCon] = fwDStepCon(X0, V);
    a = aGm + (alp - .5) * aCon;
    b = bGm + (alp - .5) * bCon;
    t = fwStepOpt(a, b);

    % update
    X = X0 + t * V;

    % debug
    if isObj
        objs(nIt) = pathDObj(X, alp);
    end

    % stop condition
    if norm(X(:) - X0(:)) < eps || t < eps
        break;
    end

    % store
    X0 = X;
    
    % debug
%     if 1
%         figure(1);
%         clf;
%         subplot(1, 2, 1);
%         cla;
%         shM(X);
%         
%         subplot(1, 2, 2);
%         plot(1 : nIt, ones(1, nIt), 'ro-');
%         drawnow;        
%     end
end
