function [asg, X0] = mysmfp(K, K1, K2, K3, ns, asgT)
% Spectral matching with integer projection.
%
% Reference
%   M. Leordeanu and M. Hebert and R. Sukthankar, "An Integer Projected
%   Fixed Point Method for Graph Matching and MAP Inference", in NIPS, 2009
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
%   K       -  affinity matrix, nn x nn (sparse)
%   ns      -  #nodes, 1 x 2
%   asgT    -  ground-truth assignment (can be [])
%
% Output
%   asg     -  assignment
%     alg   -  'smfp'
%     X     -  permutation matrix, n1 x n2
%     acc   -  accuracy
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-25-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-19-2011

% function parameter
%parIni = st('alg', 'sm', 'top', 'pow');
parIni = st('alg', 'unif', 'nor', 'n');
parPosC = st('alg', 'ipfpmy', 'deb', 'y');
parPosD = st('alg', 'gre');
prIn('mysmfp', 'ini %s, posC %s, posD %s', parIni.alg, parPosC.alg, parPosD.alg);

% initialization
X0 = asgIni(K, ns, parIni);

% continous -> discrete
XC = asgPosDIpfpMy(K1, K2, K3, X0, parPosC);

% continous -> discrete
X = asgPosD([], XC, parPosD);

% compare with ground-truth
acc = matchAsg(X, asgT);

% store
asg.alg = 'mysmfp';
asg.X = X;
asg.acc = acc;

prOut;
