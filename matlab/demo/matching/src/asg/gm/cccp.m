function [X, obj, objs] = cccp(K1, K2, X0, par)
% Convex-concave procedure.
%
% Remark
%   nn = n1 x n2
%
% Input
%   K1       -  1st affinity matrix, nn x nn (sparse)
%   K2       -  2nd affinity matrix, nn x nn (sparse)
%   X0       -  initial assignment, n1 x n2
%   par      -  parameter
%     nItMa  -  #maximum iteration steps, {50}
%
% Output
%   X        -  permutation matrix, n1 x n2
%
% History
%   create   -  Marius Leordeanu (leordeanu@gmail.com), 02-25-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-21-2011

% function parameter
nItOutMa = ps(par, 'nItOutMa', 10);
nItInMa = ps(par, 'nItInMa', 50);
isDeb = psY(par, 'deb', 'n');
nItAl = ps(par, 'nItAl', 10);

% dimension
[n1, n2] = size(X0);
ns = [n1 n2];

if isDeb
    rows = 1; cols = 3;
    AxOut = iniAx(10, rows, cols, [250 * rows, 250 * cols], 'pos', [0 0 .8 1]);
    haOut = [];
    
    rows = 1; cols = 4;
    AxIn = iniAx(11, rows, cols, [250 * rows, 250 * cols], 'pos', [0 0 .8 1]);
    haIn = [];
end

% initial
x0 = X0(:);
als = linspace(0, 1, nItAl);

% convex-concave procedure
[objOuts] = zeross(1, nItOutMa);
[objIns, lambdas] = zeross(1, nItInMa);

for nItOut = 1 : nItOutMa
    l = K1 * x0;

    % Frank-Wolfe algorithm
    xx0 = x0;
    
    for nItIn = 1 : nItInMa
        a = l - K2 * xx0;
        A = reshape(a, ns);

        % direction
        A = max(A(:)) - A;
        X = hungarian(A);
        x = X(:);
        dx = x - xx0;
        
        % step size
        C = dx' * (K2 * xx0 - l);
        D = dx' * K2 * dx;
        lambda = min([1, -C / D]);

        % update
        xx = xx0 + lambda * dx;
        
        % objective
        lambdas(nItIn) = lambda;
        objIns(nItIn) = xx' * K2 * xx - 2 * l' * xx;
        
        % debug
        if 0
            haIn = debIn(haIn, AxIn, nItIn, objIns, xx0, xx, lambdas, ns);
        end

        % stop condition
        if norm(xx - xx0) < eps
            break;
        end

        % store
        xx0 = xx;
    end
    x = xx;
    
    % objective
    objOuts(nItOut) = x' * (K1 - K2) * x;

    % debug
    if isDeb
        haOut = debOut(haOut, AxOut, nItOut, objOuts, x0, x, ns);
    end
    
    % stop condition
    if norm(x - x0) < eps
        break;
    end
    
    % store
    x0 = x;
end
obj = objOuts(nItOut);
objs = objOuts(1 : nItOut);

% reshape
X = reshape(x, [n1 n2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ha = debOut(ha, Ax, nIt, objs, x0, x, ns)

% score
shIt(objs(1 : nIt), ones(1, nIt), 'ax', Ax{1, 1}, 'mkSiz', 7, 'itMa', 0); 
title('objOut');

% initial
X0 = reshape(x0, ns);
if nIt == 1
    shM(X0, 'ax', Ax{1, 2});
end

% x
X = reshape(x, ns);
if nIt == 1
    ha.hX = shM(X, 'ax', Ax{1, 3});
else
    shMUpd(ha.hX, X);
end

drawnow;
pause(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ha = debIn(ha, Ax, nIt, objs, x0, x, lambdas, ns)

% score
shIt(objs(1 : nIt), ones(1, nIt), 'ax', Ax{1, 1}, 'mkSiz', 7, 'itMa', 0); 
title('objIn');

% initial
X0 = reshape(x0, ns);
if nIt == 1
    shM(X0, 'ax', Ax{1, 2});
end

% x
X = reshape(x, ns);
if nIt == 1
    ha.hX = shM(X, 'ax', Ax{1, 3});
else
    shMUpd(ha.hX, X);
end

% step size
shIt(lambdas(1 : nIt), ones(1, nIt), 'ax', Ax{1, 4}, 'mkSiz', 7, 'itMa', 0); 
title('step size');

drawnow;
pause(1);

