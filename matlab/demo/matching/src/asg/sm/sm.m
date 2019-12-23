function asg = sm(K, ns, asgT)
% Spectral matching for quadratic assignment.
%
% Reference
%   M. Leordeanu and M. Hebert, "A Spectral Technique
%   for Correspondence Problems Using Pairwise Constraints", in ICCV, 2005
%
% Math
%   This algorithm is to obtain the optimal P for the following problem
%     max_x   x' * K * x
%     s.t.    A * x <= 1
%
% Remark
%   nn = n1 x n2
%
% Input
%   K        -  affinity matrix, nn x nn (sparse)
%   ns       -  #nodes, 1 x 2
%   asgT     -  ground-truth assignment (can be [])
%   parIni   -  parameter for initialization
%   parPosC  -  parameter for continuous-continuous post-processing
%   parPosD  -  parameter for continuous-discrete post-processing
%
% Output
%   asg     -  assignment
%     alg   -  'sm'
%     X     -  permutation matrix, n1 x n2
%     acc   -  accuracy
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-25-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-19-2011

% function parameter
parIni = st('alg', 'sm', 'top', 'eigs');
parPosC = st('alg', 'none');
parPosD = st('alg', 'hun');
prIn('sm', 'ini %s, posC %s, posD %s', parIni.alg, parPosC.alg, parPosD.alg);

% initialization
X0 = asgIni(K, ns, parIni);

% continous -> continous
XC = asgPosC([], X0, parPosC);

% continous -> discrete
X = asgPosD([], XC, parPosD);

% compare with ground-truth
acc = matchAsg(X, asgT);

% store
asg.alg = 'sm';
asg.X = X;
asg.acc = acc;

prOut;
