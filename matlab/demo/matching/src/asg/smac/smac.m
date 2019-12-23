function asg = smac(K, ns, asgT)
% Spectral graph matching with affine constraint (SMAC).
%
% Reference
%   Timothee Cour, Praveen Srinivasan, Jianbo Shi, 
%   "Balanced Graph Matching", in NIPS, 2008
%
% Math
%   This algorithm is to obtain the optimal x for the following problem
%     max_x   x' * K * x
%     s.t.    A * x <= 1
%
% Remark
%   nn = n1 x n2
%
% Input
%   K       -  affinity matrix, nn x nn (sparse)
%   ns      -  #nodes, 1 x 2
%   asgT    -  ground-truth assignment (can be [])
%
% Output
%   asg     -  assignment
%     alg   -  'smac'
%     X     -  permutation matrix, n1 x n2
%     acc   -  accuracy
%
% History
%   create  -  Timothee Cour (timothee.cour@gmail.com), 02-01-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-15-2011

% function parameter
parIni = st('alg', 'smac');
parPosC = st('alg', 'grad');
parPosD = st('alg', 'hun');
prIn('smac', 'ini %s, posC %s, posD %s', parIni.alg, parPosC.alg, parPosD.alg);

% initialization
X0 = asgIni(K, ns, parIni);

% continous -> continous
XC = asgPosC(K, X0, parPosC);

% continous -> discrete
X = asgPosD([], XC, parPosD);

% post-processing again
% X = compute_ICM_graph_matching(X, E12, K);

% compare with ground-truth
acc = matchAsg(X, asgT);

% store
asg.alg = 'smac';
asg.X = X;
asg.acc = acc;

prOut
