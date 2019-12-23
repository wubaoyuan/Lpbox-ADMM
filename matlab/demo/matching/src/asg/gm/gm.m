function asg = gm(K, Ct, asgT, parIni, parPosC, parPosD)
% Graph matching.
%
% This function can be used as the interface of the following algorithms:
%   Graudate Assignment (GA)
%   Spectral Matching (SM)
%   Spectral Matching with Affine Constraint (SMAC)
%   Integer Projected Fixed Point (IPFP)
%   Reweighted Random Walks Matching (RRWM)
%
% Math
%   This code is to solve the following problem:
%     max_X   vec(X)' * K * vec(X)
%     s.t.    X is a permutation matrix
%
% Remark
%   nn = n1 x n2
%
% Input
%   K        -  affinity matrix, nn x nn (sparse)
%   Ct       -  correspondence constraint, n1 x n2
%                 Ct_ij = 1: i and j can be matched
%                 Ct_ij = 0: i and j cannot be matched
%   asgT     -  ground-truth assignment (can be [])
%   parIni   -  parameter for initialization
%   parPosC  -  parameter for continuous-continuous post-processing
%   parPosD  -  parameter for continuous-discrete post-processing
%
% Output
%   asg      -  assignment
%     alg    -  algorithm name
%     X      -  binary correspondence matrix, n1 x n2
%     acc    -  accuracy (= 0 if asgT is [])
%     obj    -  objective value
%     tim    -  time cost
%
% History    
%   create   -  Feng Zhou (zhfe99@gmail.com), 01-25-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 05-07-2013

% function parameter
prIn('gm', 'ini %s, posC %s, posD %s', parIni.alg, parPosC.alg, parPosD.alg);
ha = tic;

% initialization
X0 = gmIni(K, Ct, parIni);

% continous -> continous
XC = gmPosC(K, Ct, X0, parPosC);

% continous -> discrete
X = gmPosD(K, Ct, XC, parPosD);

% compare with ground-truth
acc = matchAsg(X, asgT);
%acc = matchAsg(X(1:15, 1:15), asgT);

% store
asg.alg = sprintf('gm+%s+%s+%s', parIni.alg, parPosC.alg, parPosD.alg);
asg.X = X;
asg.acc = acc;
asg.obj = X(:)' * K * X(:);
asg.tim = toc(ha);

prOut;
