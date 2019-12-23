function asg = laicClus(gph1ss, wClus, Pt2, C, Ct, par)
% An extension of LAIC algorithm for matching with several template graphs.
%
% References
%   H. Li, E. Kim, X. Huang, and L. He, 
%   "Object Matching Using a Locally Affine-Invairant Constraint", in CVPR, 2010.
%
% Input
%   gph1s   -  template graph, 1 x m (cell)
%   Pt2     -  testing point set, 2 x n
%   C       -  cost matrix, k x n
%   Ct      -  constraint, k x n
%   par     -  parameter
%     lam   -  The parameter weighting feature cost and geometric cost, {1}
%
% Output
%   asg
%     alg   -  'laic'
%     X     -  correspondence Hmatrix, k x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-02-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 08-26-2012

% function parameters
lam = ps(par, 'lam', 10);
rat = ps(par, 'rat', .7);
nItMa = ps(par, 'nItMa', 0);
isDeb = psY(par, 'deb', 'n');
prIn('laicClus', 'lam %.2f, rat %.2f, nItMa %d', lam, rat, nItMa);

% dimension
[k, n] = size(C);
kClu = length(gph1ss);

% weights
weis = par.weis;

% debug
if isDeb
    Rs = par.Rs;
    PtT = par.PtT;
    rows = 4; cols = 8;
    Ax = iniAx(10, rows, cols, [200 * rows, 200 * cols]);
    ha = [];
end

% template graph
Hs = cell(1, kClu);
for cClu = 1 : kClu
    gph1s = gph1ss{cClu};
    m = length(gph1s);
    
    if m == 0
        continue;
    end

    [Pts, As] = cellss(1, m);
    for i = 1 : m
        [Pt1, Eg1] = stFld(gph1s{i}, 'Pt', 'Eg');

        % edge -> adjacency matrix
        A0 = gphEg2Adj(Eg1, k);
        
        % make sure every node has at least three neighbours
        A = laicValA(A0);
        
        Pts{i} = Pt1;
        As{i} = A;
    end

    % calculate the reconstruction matrix
    Hs{cClu} = calcReconCoes(Pts, As);
end

% basis & region
[ba0s, Reg0] = laicBas(Pt2, Ct);

% LP relaxation
XC0 = laicLPClus(Hs, wClus, Pt2', C, ba0s, lam, weis);

% ICM
XD0 = laicICMClus(XC0, Hs, wClus, Pt2', C, ba0s, lam, weis);

if isDeb
    ha = deb(ha, Ax, Rs, Pt2, ba0s, Reg0, XC0, XD0, PtT);
end
XC = XC0;
XD = XD0;

% shrink the trust-region
prCIn('shrink', nItMa, 1);
for nIt = 1 : nItMa
    prC(nIt);
    
    % update region
    [bas, Reg] = laicRegNew(X, Pt2, ba0s, Reg0, rat);
    
    % compute the LP relaxation
    XC = laicLPs(Hs, Pt2', C, bas, lam, weis);
    
    % ICM
    XD = laicICM(XC, Hs, Pt2', C, bas, lam, weis);
    
    if isDeb
        ha = deb(ha, Ax, nIt, F, X, Pt2, bas, Reg);
    end

    % store
    Reg0 = Reg;
end
prCOut(nIt);

% obj
[obj, objs, obj1s] = laicObjClus(XD, Pt2, C, Hs, wClus, lam, weis);

% sort point
[~, idxTrs] = sort(obj1s);

% store
asg.alg = 'laic2s';
asg.X = XD;
asg.XC = XC;
asg.PtNew = Pt2 * XD';
asg.obj = obj;
asg.objs = objs;
asg.C = C;
asg.Hs = Hs;
asg.Pt2 = Pt2;
asg.idxTrs = idxTrs;
asg.lam = lam;
asg.weis = weis;

prOut;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ha = deb(ha, Ax, Rs, Pt, bas, Reg, XC, XD, PtT)
% Debug

% label name
nms = lfpwLMNmIdx(2);

% dimension
k = length(bas);
[rows, cols] = size(Ax);

% per region
for c = 1 : k
    Rc = Rs(:, :, c);
    idx = bas{c};
    
    Ptc = Pt(:, idx);
    xC = XC(c, idx);
    xD = XD(c, idx);

    mi = min([Ptc, PtT(:, c)], [], 2);
    ma = max([Ptc, PtT(:, c)], [], 2);

    % show image
    shM(Rc, 'ax', Ax{c});
    plot(Ptc(1, :), Ptc(2, :), 'ro', 'markersize', 7, 'linewidth', 1, 'markeredgecolor', 'k');
    set(gca, 'xlim', [mi(1) ma(1)], 'ylim', [mi(2) ma(2)]);
    
    % plot region
    miX = Reg(1, 1, c);
    miY = Reg(2, 1, c);
    maX = Reg(1, 2, c);
    maY = Reg(2, 2, c);
    plot([miX miX maX maX miX], [miY maY maY miY miY], 'b-');
    
    % ground-truth
    plot(PtT(1, c), PtT(2, c), 'rx', 'linewidth', 1, 'markersize', 7);

    % continuous solution
    idxC = find(xC > 0);
    PtC = Ptc * xC';
    plot(PtC(1), PtC(2), 'b+', 'linewidth', 1, 'markersize', 7);    
    plot(Ptc(1, idxC), Ptc(2, idxC), 'b+', 'markersize', 7);
    
    % discrete solution
    idxD = find(xD > 0);
    PtD = Ptc * xD';
    plot(PtD(1), PtD(2), 'r+', 'linewidth', 1);    

    title(nms{c});    
end
for c = k + 1 : rows * cols
    set(Ax{c}, 'visible', 'off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bas, Reg] = laicBas(Pt, Ct)

[k, n] = size(Ct);

bas = cellss(k, 1);
Reg = zeros(2, 2, k);
for c = 1 : k
    bas{c} = find(Ct(c, :) == 1);
    
    Ptc = Pt(:, bas{c});

    boxMi = min(Ptc, [], 2);
    boxMa = max(Ptc, [], 2);
    
    Reg(:, :, c) = [boxMi, boxMa];
end
