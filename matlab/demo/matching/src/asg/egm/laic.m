function asg = laic(gph1, Pt2, C, XT, par)
% LAIC algorithm.
%
% References
%   H. Li, E. Kim, X. Huang, and L. He, 
%   "Object Matching Using a Locally Affine-Invairant Constraint", in CVPR, 2010.
%
% Input
%   gph1    -  template graph
%   Pt2     -  testing point set, 2 x n
%   C       -  cost matrix, k x n
%   XT      -  ground-truth correspondence, k x n | []
%   par     -  parameter
%     lam   -  The parameter weighting feature cost and geometric cost, {1}
%
% Output
%   asg
%     alg   -  'laic'
%     X     -  correspondence Hmatrix, k x n
%
% History
%   create  -  Hongsheng Li (h.li@lehigh.edu), 12-13-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-27-2012

% function parameters
lam = ps(par, 'lam', 1);
prIn('laic', 'lam %.2f', lam);

% dimension
[k, n] = size(C);

% template graph
[Pt1, Eg1] = stFld(gph1, 'Pt', 'Eg');

% edge -> adjacency matrix
A0 = gphEg2Adj(Eg1, k);

% make sure every node has three neighbours
A = laicValA(A0);

% calculate the reconstruction matrix
H = calcReconCoe(Pt1', A);

% basis
bas = laicHull(Pt2', C);

% LP relaxation
[X, objTmp, res1, res2, Aeq, beq] = laicLP(H, Pt2', C, bas, lam);

% test
% [res1a, res2a, XPosa, XNega] = laicLPX2Res(XT, H, Pt2', bas);
% beqT = Aeq * [res1a; res2a];
% equal('b', beqT, beq);
% equal('XPos + XNeg', XPosa - XNega, H * XT * Pt2');

% objective
[objT, objT1, objT2] = laicObj(XT, Pt2, C, H, lam);
[obj, obj1, obj2, obj1s] = laicObj(X, Pt2, C, H, lam);

% sort the confidence of node
[~, idxTrs] = sort(obj1s);

% store
asg.alg = 'laic';
asg.X = X;
asg.idxTrs = idxTrs;
asg.objT = objT;
asg.objT1 = objT1;
asg.objT2 = objT2;
asg.obj = obj;
asg.obj1 = obj1;
asg.obj2 = obj2;
asg.PtNew = Pt2 * X';

prOut;
