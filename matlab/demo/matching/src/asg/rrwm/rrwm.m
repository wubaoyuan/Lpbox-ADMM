function asg = rrwm(K, ns, asgT)
% Reweighted Random Walk Matching.
%
% Reference
%   M. Cho, J. Lee, and K. M. Lee. "Reweighted random walks for
%   graph matching", In ECCV, 2010. 
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
%   n0      -  #destined pairs
%   asgT    -  ground-truth assignment (can be [])
%
% Output
%   asg     -  assignment
%     alg   -  'rrwm'
%     X     -  permutation matrix, n1 x n2
%     acc   -  accuracy
%
% History
%   create  -  Minsu Cho (chominsu@gmail.com), 09-25-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-15-2011

% function parameter
parIni = st('alg', 'unif', 'nor', 'n');
parPosC = st('alg', 'rrwm');
parPosD = st('alg', 'gre');
prIn('rrwm', 'ini %s, posC %s, posD %s', parIni.alg, parPosC.alg, parPosD.alg);

% initialization
X0 = asgIni([], ns, parIni);

% continous -> continous
XC = asgPosC(K, X0, parPosC);

% continous -> discrete
X = asgPosD([], XC, parPosD);

% compare with ground-truth
acc = matchAsg(X, asgT);

% store
asg.alg = 'rrwm';
asg.X = X;
asg.acc = acc;

prOut;
