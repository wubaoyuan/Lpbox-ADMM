function asg = ga(K, ns, asgT)
% Graduate assignment for graph matching.
%
% Reference
%   S. Gold and A. Rangarajan. "A graduated assignment algorithm
%   for graph matching", IEEE Trans. Pattern Anal. Mach. Intell.,
%   18(4):377–388, 1996.
%
% Math
%   This algorithm is to obtain the optimal X for the following problem
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
%     alg   -  'ga'
%     X     -  permutation matrix, n1 x n2
%     acc   -  accuracy
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 10-01-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-19-2011

% function parameter
parIni = st('alg', 'unif', 'nor', 'n');
parPosC = st('alg', 'grad', 'b0', .5, 'bMax', 10);
% parPosC = st('alg', 'grad', 'b0', max(ns), 'bMax', 200); % better
parPosD = st('alg', 'gre');
prIn('ga', 'ini %s, posC %s, posD %s', parIni.alg, parPosC.alg, parPosD.alg);

% initialization
X0 = asgIni([], ns, parIni);

% continous -> continous
XC = asgPosC(K, X0, parPosC);

% continous -> discrete
X = asgPosD([], XC, parPosD);

% compare with ground-truth
acc = matchAsg(X, asgT);

% store
asg.alg = 'ga';
asg.X = X;
asg.acc = acc;

prOut
