function [X, X2, obj, nIts, Xs, objs, objGms, objCons, objVexs, objCavs, useIps, objInss, objIn2ss] ...
    = pathDIter(alps, nItMa, nHst, isIp, isDeb, isDis, Y, isYM)
% Path-following iteration for optimizing the factorized graph matching.
%
% Remark
%   The edge is directed and the edge feature is asymmetric. 
%
% Input
%   alps      -  weights, 1 x nAlp
%   nItMa     -  #maximum iteration steps for each scale of alpha
%   nHst      -  #saved steps in the modified Frank-Wolfe (MFW) iteration
%   isIp      -  flag of using IPFP to adjust the searching direction, 0 | 1
%   isDeb     -  flag of debugging, 0 | 1
%   isDis     -  flag of discretizing the result after the last step, 0 | 1
%   Y         -  specified initial correspondence, n1 x n2 | []
%   isYM      -  flag of starting from the middle, 0 | 1
%             
% Output      
%   X         -  correspondence matrix, n1 x n2
%   X2        -  correspondence matrix, n1 x n2
%   obj       -  objective
%
%   ** For debugging **
%   nIts      -  #iterations, 1 x nAlp
%   Xs        -  correspondence matrices, 1 x nAlp (cell), n1 x n2
%   objs      -  objective at each iteration, 1 x nAlp
%   objGms    -  objective (J_gm), 1 x nAlp
%   objCons   -  objective (J_con), 1 x nAlp
%   objVexs   -  objective (J_vex), 1 x nAlp
%   objCavs   -  objective (J_cav), 1 x nAlp
%   useIps    -  state of using IPFP, 1 x nAlp
%   objInss   -  objecitve at each step, 1 x nAlp (cell)
%   objIn2ss  -  objecitve at each step, 1 x nAlp (cell)
%             
% History     
%   create    -  Feng Zhou (zhfe99@gmail.com), 09-01-2011
%   modify    -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% global (specified in the pathDStore.m)
global Ct ns XT;

% initialize from a doubly stochastic matrix
X0 = gmIniUnif(Ct, st('nor', 'doub'));

% starting from the middle with the specified initialization Y
if ~isempty(Y) && isYM
    head = find(alps == .5);
    alps = alps(head : end);
end
objGm0 = pathDObjGm(X0);

% dimension
nAlp = length(alps);
[nIts, objs, objGms, objCons, objVexs, objCavs, useIps] = zeross(1, nAlp);
[Xs, objInss, objIn2ss] = cellss(1, nAlp);

% figure for debugging
if isDeb
    rows = 1; cols = 4;
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
    [X, nIts(iAlp), objIns] = mfwDIter(X0, alp, nItMa, nHst, isDeb);    

    % compared to FW (only for debugging)
    if isDeb
        [~, ~, objIn2s] = fwDIter(X0, alp, nItMa, isDeb);
    else
        objIn2s = objIns;
    end

    % objective
    [objs(iAlp), objGms(iAlp), objCons(iAlp), objVexs(iAlp), objCavs(iAlp)] = pathDObj(X, alp);
    
    % using IPFP to refine the result
    if isIp && ((iAlp > 1 && objGms(iAlp) < objGms(iAlp - 1)) || (iAlp == 1 && objGm0 > objGms(iAlp)))
        X = ipfpAIter(X0, 1);
        useIps(iAlp) = 1;
        pr('using ipfp');
    end

    % refine: using specified initial correspondence
    if iAlp < nAlp && alps(iAlp + 1) == .5 && ~isempty(Y)
        % objective
        objGmX = pathDObjGm(X);
        objGmY = pathDObjGm(Y);

        if objGmX < objGmY
            X = Y;
            pr('using Y');
        end
    end

    % compute objective for ground-truth
    if isempty(XT)
        objT = [];
        objGmT = [];
    else
        [objT, objGmT] = pathDObj(XT, alp);
    end

    % save X and obj
    if isDeb
        Xs{iAlp} = X;
        objInss{iAlp} = objIns;
        objIn2ss{iAlp} = objIn2s;
    end

    % debug
    if isDeb
        ha = deb(ha, Ax, iAlp, X, objs, objGms, objIns, objIn2s, objT, objGmT, useIps);
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
% obj = pathDObj(X, 1);
[~, obj] = pathDObj(X, 1);

% re-size to the original size
X2 = X;
X = X(1 : ns(1), 1 : ns(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ha = deb(ha, Ax, iAlp, X, objs, objGms, objIns, objIn2s, objT, objGmT, useIps)
% Debug.

% correpondence
shM(X, 'ax', Ax{1, 1});
title('X');

shIt(objs(1 : iAlp), ones(1, iAlp), 'ax', Ax{1, 2}, 'itMa', 0);
title('obj (path)');
mi = min(objs(1 : iAlp));
ma = max(objs(1 : iAlp));
for i = 1 : iAlp
    if useIps(i)
        plot([i i], [mi ma], '--k');
    end
end

if ~isempty(objT)
    plot([1, length(objs)], [0 0] + objT, '-r');
    objss = [objs, objT];
    if abs(min(objss) - max(objss)) > eps
        set(gca, 'ylim', [min(objss), max(objss)]);
    end
end

shIt(objGms(1 : iAlp), ones(1, iAlp), 'ax', Ax{1, 3}, 'itMa', 0);
if ~isempty(objGmT)
    plot([1, length(objs)], [0 0] + objGmT, '-r');
    objGmss = [objGms, objGmT];
    if abs(min(objGmss) - max(objGmss)) > eps
        set(gca, 'ylim', [min(objGmss), max(objGmss)]);    
    end
end
title('obj (gm)');

set(gcf, 'CurrentAxes', Ax{1, 4});
cla;
hold on;
plot(objIn2s, '-sb');
plot(objIns, '-or');
legend('FW', 'MFW');

drawnow;
