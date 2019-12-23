function [X, obj, nIts, Xs, objs, objGms, objCons, objVexs, objCavs] ...
    = pathUIter(ns, alps, nItMa, isObj, isSvX, isDeb, isDis, isIp, nHst)
% Path-following iteration for optimizing the factorized graph matching
% with an asymmetric edge feature.
%
% Input
%   ns       -  #nodes, 1 x 2
%   alps     -  weights, 1 x nAlp
%   nItMa    -  #maximum iteration steps for each scale of alpha
%   isObj    -  flag of computing objective, 0 | 1
%   isSvX    -  flag of saving X, 0 | 1
%   isDeb    -  flag of debugging, 0 | 1
%   isDis    -  flag of discreting the result, 0 | 1
%   isIp     -  flag of using IPFP to refine the result, 0 | 1
%
% Output
%   X        -  correspondence matrix, n1 x n2
%   obj      -  objective
%   nIts     -  #iterations, 1 x nAlp
%   Xs       -  correspondence matrices, 1 x nAlp (cell), n1 x n2
%   objs     -  objective, 1 x nAlp
%   objGms   -  objective, 1 x nAlp
%   objCons  -  objective, 1 x nAlp
%   objVexs  -  objective, 1 x nAlp
%   objCavs  -  objective, 1 x nAlp
%
% History   
%   create   -  Feng Zhou (zhfe99@gmail.com), 09-01-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% dimension
n = max(ns);

% global variables;
global Ct;

% initialize
X0 = gmIniUnif(Ct, st('nor', 'doub'));
X00 = X0;

% dimension
nAlp = length(alps);

% allocate space for saving variables
[nIts, objs, objGms, objVexs, objCavs] = zeross(1, nAlp);
Xs = cell(1, nAlp);

% figure for debugging
if isDeb
    rows = 1; cols = 2;
    Ax = iniAx(11, rows, cols, [300 * rows, 300 * cols], 'pos', [0 0 1 1]);
    ha = [];
end

% path-following
prCIn('path', nAlp, .1);
for iAlp = 1 : nAlp
    prC(iAlp);

    % scale of alpha
    alp = alps(iAlp);

    % MFW
    [X, nIts(iAlp)] = mfwUIter(X0, alp, nItMa, nHst, isObj);
    
    % IPFP
    if isIp
        % objective
        [objs(iAlp), objGms(iAlp)] = pathUObj(X, alp);

        % using IPFP to refine the result
        if iAlp > 1 && objGms(iAlp) < objGms(iAlp - 1)
            X = ipfpSIter(X0, 1);
        end
    end

    % save objective
    if isObj || isDeb
        [objs(iAlp), objGms(iAlp), objVexs(iAlp), objCavs(iAlp)] = pathUObj(X, alp);
    end

    % save X
    if isSvX
        Xs{iAlp} = X;
    end

    % debug
    if isDeb
        ha = deb(ha, Ax, iAlp, objs, X);
    end

    % store
    X0 = X;
end
prCOut(nAlp);
pr('FW iter: %d', sum(nIts));

% post-processing
if isDis
    XC = X;
    X = gmPosDHun(X);
    if ~equal('XC', XC, X, 'pr', 'n')
        pr('non-discrete');
    end
end

% objective
obj = pathUObj(X, 1);

% re-size to the original size
X = X(1 : ns(1), 1 : ns(2));

Xs{1} = X00;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ha = deb(ha, Ax, iAlp, objs, X)
% Debug.

shIt(objs(1 : iAlp), ones(1, iAlp), 'ax', Ax{1, 1}, 'itMa', 0);

% correpondence
shM(X, 'ax', Ax{1, 2});

drawnow;
