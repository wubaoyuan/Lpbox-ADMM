function asg = icp(gphs, tran0, asgT, parIcp, parTran)
% Iterative cloest point (ICP).
%
% Input
%   gphs       -  graphs, 1 x 2 (cell)
%     Pt       -  node position, d x ni
%   tran0      -  initial transformation
%   asgT       -  ground-truth assignment (can be [])
%   parIcp     -  parameter for ICP
%     nItEmMa  -  #maximum EM iteration steps, {100}
%     deb      -  flag of debugging, 'y' | {'n'}
%     th       -  threshold for stop, {1e-6}
%   parTran    -  parameter for transformation
%     algT     -  transformation type, {'sim'} | 'aff' | 'non'
%     sigW     -  sigma for computing RBF kernel
%     P        -  basis point for computing RBF kernel
%              
% Output       
%   asg        -  assignment
%     alg      -  'icp'
%     X        -  correspondence matrix, n1 x n2
%     tran     -  transformation
%     acc      -  accuracy
%     obj      -  objective
%
% History
%   create     -  Feng Zhou (zhfe99@gmail.com), 02-16-2012
%   modify     -  Feng Zhou (zhfe99@gmail.com), 12-11-2012

% function parameter
algT = ps(parTran, 'algT', 'sim');
nItEmMa = ps(parIcp, 'nItEmMa', 100);
isDeb = psY(parIcp, 'deb', 'n');
th = ps(parIcp, 'th', 1e-6);
prIn('icp', 'algT %s, nItEmMa %d', algT, nItEmMa);

% dimension
P1 = gphs{1}.Pt;
P2 = gphs{2}.Pt;
n1 = size(P1, 2);
n2 = size(P2, 2);
gph0s = gphs;

% figure for debugging
if isDeb
    rows = 2; cols = 3;
    rows = 1; cols = 4;    
    Ax = iniAx(10, rows, cols, [300 * rows, 300 * cols], 'pos', [0 0 1 1]);
    ha = [];
    set(Ax{1, 3}, 'visible', 'off');
end

% EM iteration
[objs, its] = zeross(1, nItEmMa);
prCIn('EM', nItEmMa, 1);
for nItEm = 1 : 2 : nItEmMa
    prC(nItEm);

    % optimal correspondence
    X = corrOpt(P1, P2, tran0);
    its(nItEm) = 1;
    objs(nItEm) = icpObj(gphs, X, tran0);

    % debug
    if isDeb
        ha = deb(ha, Ax, nItEm, objs, its, gphs, gph0s, X, tran0);
    end

    % optimal transformation
    tran = tranOpt(gphs, 0, X, parTran);
    its(nItEm + 1) = 2;
    objs(nItEm + 1) = icpObj(gphs, X, tran);

    % debug
    if isDeb
        ha = deb(ha, Ax, nItEm + 1, objs, its, gphs, gph0s, X, tran);
    end

    % stop condition
    if tranDif(tran0, tran, st('th', th))
        break;
    end

    % store
    tran0 = tran;
end
prCOut(nItEm + 1);
objs(nItEm + 2 : end) = objs(nItEm + 1);

% obj
obj = icpObj(gphs, X, tran);

% matching with ground-truth assignment if possible
acc = matchAsg(X, asgT);

% apply transformation
gphs{2}.Pt = tranRun(gphs{2}.Pt, tran);
gphs{1}.Eg = [];
gphs{2}.Eg = [];

% store
asg.alg = 'icp';
asg.X = X;
asg.tran = tran;
asg.obj = obj;
asg.objs = objs;
asg.acc = acc;
asg.gphs = gphs;

prOut;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = corrOpt(P1, P2, tran)
% Compute the optimal correspondence.
%
% Input
%   P1    -  1st point set, d x n1
%   P2    -  2nd point set, d x n2
%   tran  -  initial transformation
%         
% Output  
%   X     -  correspondence matrix, n1 x n2

% dimension
n1 = size(P1, 2);
n2 = size(P2, 2);

% transform the points
P2a = tranRun(P2, tran);

% KP
KP = 2 * P1' * P2a - (P1 .* P1)' * ones(2, n2) - ones(n1, 2) * (P2a .* P2a);

% hungrian
X = gmPosDHun(KP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ha = deb_old(ha, Ax, nItEm, objEms, itEms, gphs, gph0s, X, tran)
% Debug.

% objective
shIt(objEms(1 : nItEm), itEms(1 : nItEm), 'ax', Ax{1, 1}, 'itMa', 0);
title('objective');

% correpondence
shM(X, 'ax', Ax{1, 2});
title('correspondence matrix');

% remove edge
gphs{1}.Eg = [];
gphs{2}.Eg = [];
Pt20 = gphs{2}.Pt;

% 2nd graph after transformation
parMks = {st('cl', 'b'), st('cl', 'r')};
gphs{2}.Pt = tranRun(gphs{2}.Pt, tran);
shGph(gphs, 'ax', Ax{2, 2}, 'parMks', parMks);
axis ij;
title('current graphs');

% connection
Pt1 = gphs{1}.Pt;
Pt2 = gphs{2}.Pt;
idx = find(X);
[is, js] = ind2sub(size(X), idx);
nC = length(idx);
xs = [Pt1(1, is); Pt2(1, js); nan(1, nC)];
ys = [Pt1(2, is); Pt2(2, js); nan(1, nC)];
plot(xs(:), ys(:), '--k');
axis equal;

% original graph
parMks = {st('cl', 'b', 'lnWid', 0), st('cl', 'r', 'lnWid', 0)};
shGph(gph0s, 'ax', Ax{2, 1}, 'parMks', parMks);
axis ij;

title('original graphs');

% new graph
shGph(gphs, 'ax', Ax{2, 3}, 'parMks', parMks);
hold on;
axis ij;

[Gr2, controls] = ctps_plot_grid_gen(Pt20', 7, 5);
Gr2 = Gr2';
Gr1 = tranRun(Gr2, tran);
ctps_plot_gridbox(1, Gr1', controls, 'r', ':');
ctps_plot_gridbox(1, Gr2', controls, 'b', ':');

%set(gca, 'xlim', [mi2s(1) ma2s(1)], 'ylim', [mi2s(2) ma2s(2)]);
drawnow;
pause(1);
title('current graphs with warping grid');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ha = deb(ha, Ax, nItEm, objEms, itEms, gphs, gph0s, X, tran)
% Debug for Computational Photograph class.

% objective
shIt(objEms(1 : nItEm), itEms(1 : nItEm), 'ax', Ax{1, 1}, 'itMa', 0);
%title('objective');

% correpondence
% shM(X, 'ax', Ax{1, 2});
% title('correspondence matrix');

% remove edge
gphs{1}.Eg = [];
gphs{2}.Eg = [];
Pt20 = gphs{2}.Pt;

% 2nd graph after transformation
parMks = {st('cl', 'b'), st('cl', 'r')};
gphs{2}.Pt = tranRun(gphs{2}.Pt, tran);
shGph(gphs, 'ax', Ax{1, 3}, 'parMks', parMks);
axis ij;
%title('current graphs');

% connection
Pt1 = gphs{1}.Pt;
Pt2 = gphs{2}.Pt;
idx = find(X);
[is, js] = ind2sub(size(X), idx);
nC = length(idx);
xs = [Pt1(1, is); Pt2(1, js); nan(1, nC)];
ys = [Pt1(2, is); Pt2(2, js); nan(1, nC)];
plot(xs(:), ys(:), '--k');
axis equal;

% original graph
parMks = {st('cl', 'b', 'lnWid', 0), st('cl', 'r', 'lnWid', 0)};
shGph(gph0s, 'ax', Ax{1, 2}, 'parMks', parMks);
axis ij;
axis equal;
axis off;

%title('original graphs');

% new graph
shGph(gphs, 'ax', Ax{1, 4}, 'parMks', parMks);
hold on;
axis ij;
axis equal;
axis off;

[Gr2, controls] = ctps_plot_grid_gen(Pt20', 7, 5);
Gr2 = Gr2';
Gr1 = tranRun(Gr2, tran);
ctps_plot_gridbox(1, Gr1', controls, 'r', ':');
ctps_plot_gridbox(1, Gr2', controls, 'b', ':');

%set(gca, 'xlim', [mi2s(1) ma2s(1)], 'ylim', [mi2s(2) ma2s(2)]);
drawnow;
pause(1);
%title('current graphs with warping grid');
