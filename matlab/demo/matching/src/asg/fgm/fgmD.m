function asg = fgmD(KP, KQ, Ct, gphs, asgT, par)
% Factorized graph matching.
%
% Remark
%   The edge is directed and the edge feature is asymmetric.
%
% Reference
%   F. Zhou and F. De la Torre, "Deformable Graph Matching", in CVPR, 2013.
%
% Input
%   KP       -  node affinity matrix, n1 x n2
%   KQ       -  edge affinity matrix, m1 x m2
%   Ct       -  correspondence constraint, n1 x n2
%                 Ct_ij = 1: i and j can be matched
%                 Ct_ij = 0: i and j cannot be matched
%   gphs     -  graphs, 1 x 2 (cell)
%     G      -  node-edge adjacency (for starting point), n x mi
%     H      -  node-edge adjacency (for ending point), n x mi
%   par      -  parameter
%     nAlp   -  #alpha, {100}
%     nItMa  -  #maximum iteration steps for each scale of alpha, {100}
%     nHst   -  #history nodes for modifed FW algorithm, {10}
%     ip     -  flag of using IPFP to improve the algorithm, {'y'} | 'n'
%     deb    -  flag of debugging, 'y' | {'n'}
%
% Output
%   asg      -  assignment
%     alg    -  'fgmD'
%     X      -  correspondence matrix, n1 x n2
%     acc    -  accuracy
%     obj    -  objective
%
% History    
%   create   -  Feng Zhou (zhfe99@gmail.com), 09-01-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% function parameter
nAlp = ps(par, 'nAlp', 101);
nItMa = ps(par, 'nItMa', 100);
nHst = ps(par, 'nHst', 10);
isIp = psY(par, 'ip', 'n');
isDeb = psY(par, 'deb', 'n');
prIn('fgmD', 'nAlp %d, nItMa %d, nHst %d, isIp %d', ...
     nAlp, nItMa, nHst, isIp);

% weight
alps = linspace(0, 1, nAlp);

% store variables
XQ1 = ps(par, 'XQ1', []);
XQ2 = ps(par, 'XQ2', []);
pathDStore(KP, KQ, Ct, XQ1, XQ2, gphs, 'path', []);

% path-following
[X, ~, obj, nIts, Xs, objs, objGms, objCons, objVexs, objCavs, useIps, objInss, objIn2ss] = ...
    pathDIter(alps, nItMa, nHst, isIp, isDeb, 1, [], 0);

% matching with ground-truth assignment if possible
acc = matchAsg(X, asgT);

% store
asg.alg = 'fgmD';
asg.X = X;
asg.obj = full(obj);
asg.acc = acc;

% for debug
asg.nIts = nIts;
asg.Xs = Xs;
asg.objs = objs;
asg.objGms = objGms;
asg.objVexs = objVexs;
asg.objCavs = objCavs;
asg.objCons = objCons;
asg.objInss = objInss;
asg.objIn2ss = objIn2ss;
asg.useIps = useIps;

prOut;
