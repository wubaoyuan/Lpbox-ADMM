function asg = laics(Hs, Pt, C, Ct, lams, par)
% An extension of LAIC algorithm for matching with several template graphs.
%
% References
%   H. Li, E. Kim, X. Huang, and L. He, 
%   "Object Matching Using a Locally Affine-Invairant Constraint", in CVPR, 2010.
%
% Input
%   Hs      -  reconstruct matrix, 1 x m (cell), k x k
%   Pt      -  testing point set, 2 x n
%   C       -  cost matrix, k x n
%   Ct      -  constraint, k x n
%   lams    -  weights feature cost and geometric cost, {1}
%   par     -  parameter
%     rat   -  ratio
%     nItMa -  #maximum iterations
%     deb   -  debug flag, 'y' | {'n'}
%
% Output
%   asg
%     alg   -  'laics'
%     X     -  correspondence matrix, k x n
%     XC    -  continuous correspondence, k x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-16-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 08-28-2012

% function parameters
rat = ps(par, 'rat', .7);
nItMa = ps(par, 'nItMa', 0);
isDeb = psY(par, 'deb', 'n');
prIn('laics', 'rat %.2f, nItMa %d', rat, nItMa);

% dimension
[k, n] = size(C);

% debug
if isDeb
    Rs = par.Rs;
    PtT = ps(par, 'PtT', []);
    rows = 4; cols = 8;
    Ax = iniAx(10, rows, cols, [200 * rows, 200 * cols]);
    ha = [];
end

% basis & region
[ba0s, Reg0] = laicBas(Pt, Ct);

% LP relaxation
hT = tic;
XC0 = laicLPs(Hs, Pt', C, ba0s, lams);
tiC = toc(hT);

% ICM
hT = tic;
XD0 = laicICM(XC0, Hs, Pt', C, ba0s, lams);
tiD = toc(hT);

if isDeb
    ha = deb(ha, Ax, Rs, Pt2, ba0s, Reg0, XC0, XD0, PtT);
end
XC = XC0;
XD = XD0;
bas = ba0s;

% shrink the trust-region
prCIn('shrink', nItMa, 1);
for nIt = 1 : nItMa
    prC(nIt);

    % update region
    [bas, Reg] = laicRegNew(XD0, Pt, ba0s, Reg0, rat);
    
    % compute the LP relaxation
    XC = laicLPs(Hs, Pt', C, bas, lams);
    
    % ICM
    XD = laicICM(XC, Hs, Pt', C, bas, lams);
    
    if isDeb
        ha = deb(ha, Ax, nIt, F, X, Pt, bas, Reg);
    end

    % store
    Reg0 = Reg;
    XC0 = XC;
    XD0 = XD;
end
prCOut(nIt);

% obj
[obj, objs, obj1s] = laicObjs(XD, Hs, Pt, C, bas, lams);

% sort point
[~, idxTrs] = sort(obj1s);
idxTes = zeros(1, k);
for c = 1 : k
    idxTes(c) = find(XD(idxTrs(c), :));
end

% store
asg.alg = 'laics';
asg.X = XD;
asg.XC = XC;
asg.PtNew = Pt * XD';
asg.obj = obj;
asg.objs = objs;
asg.C = C;
asg.Hs = Hs;
asg.bas = bas;
asg.Pt = Pt;
asg.idxTrs = idxTrs;
asg.idxTes = idxTes;
asg.lams = lams;
asg.tiC = tiC;
asg.tiD = tiD;

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

    mi = min(Ptc, [], 2);
    ma = max(Ptc, [], 2);
    if ~isempty(PtT)
        mi = min([mi, PtT(:, c)], [], 2);
        ma = max([ma, PtT(:, c)], [], 2);
    end

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
    if ~isempty(PtT)
        plot(PtT(1, c), PtT(2, c), 'rx', 'linewidth', 1, 'markersize', 7);
    end

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
% Obtain basis and region.

% dimension
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
