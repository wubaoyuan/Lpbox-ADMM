function asg = icpVdo(gphs, tran0, asgT, parIcp, parTran)
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
%   modify     -  Feng Zhou (zhfe99@gmail.com), 05-09-2012

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

% figure for debugging
if isDeb
    rows = 1; cols = 1;
    figSiz = [300 * rows, 300 * cols];
    Ax = iniAx(10, rows, cols, figSiz);
    ha = [];

    % output video
    path = 'icp_opt.avi';
    hw = vdoWIn(path, 'fps', 1, 'siz', figSiz, 'comp', 'vdo');
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
        ha = deb(ha, Ax, nItEm, objs, its, gphs, X, tran0);
    end

    % optimal transformation
    tran = tranOpt(gphs, 0, X, parTran);
    its(nItEm + 1) = 2;
    objs(nItEm + 1) = icpObj(gphs, X, tran);

    % debug
    if isDeb
        ha = deb(ha, Ax, nItEm + 1, objs, its, gphs, X, tran);
        hw = vdoW(hw);        
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

if isDeb
    vdoWOut(hw);
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ha = deb(ha, Ax, nItEm, objEms, itEms, gphs, X, tran)
% Debug.

gphs{1}.Eg = [];
gphs{2}.Eg = [];
parMks = {st('cl', 'b'), st('cl', 'r')};

% 2nd graph after transformation
gphs{2}.Pt = tranRun(gphs{2}.Pt, tran);
shGph(gphs, 'ax', Ax{1, 1}, 'parMks', parMks);

drawnow;
pause(1);
