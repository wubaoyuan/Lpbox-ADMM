function asg = ipfpU(KP, KQ, Ct, gphs, asgT, par)
% IPFP.
%
% Remark
%   The edge feature is symmetric.
%
% Reference
%   M. Leordeanu and M. Hebert and R. Sukthankar, "An Integer Projected
%   Fixed Point Method for Graph Matching and MAP Inference", in NIPS, 2009
%
% Input
%   KP        -  node affinity matrix, n1 x n2
%   KQ        -  edge affinity matrix, m1 x m2
%   Ct        -  constraint, n1 x n2
%                  Ct_ij = 1: i and j can be matched
%                  Ct_ij = 0: i and j cannot be matched
%   gphs      -  graphs, 1 x 2 (cell)
%     G       -  node-edge adjacency matrix, ni x mi
%     H       -  augment node-edge adjacency matrix, ni x (mi + ni)
%   asgT      -  ground-truth assignment (can be [])
%   par       -  parameter
%     nItMa   -  #maximum iteration steps for each scale of alpha, {100}
%
% Output
%   asg       -  assignment
%     alg     -  'ipfpS'
%     X       -  correspondence matrix, n1 x n2
%     acc     -  accuracy
%     obj     -  objective
%
% History     
%   create    -  Feng Zhou (zhfe99@gmail.com), 09-01-2011
%   modify    -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% function parameter
nItMa = ps(par, 'nItMa', 100);
prIn('ipfpU', 'nItMa %d', nItMa);

% dimension
n1 = size(gphs{1}.G, 1);
n2 = size(gphs{2}.G, 1);
ns = [n1, n2];

% store variables
ipfpUStore(KP, KQ, Ct, gphs);

nItMa = 20;
nHst = 5;

% dimension
n = max(ns);

% initialize
X0 = gmIniUnif(Ct, st('nor', 'doub'));

% IPFP
[X, nIt] = ipfpUIter(X0, nItMa);

% matching with ground-truth assignment if possible
acc = matchAsg(X, asgT);

% store
asg.alg = 'ipfpS';
asg.X = X;
asg.nIt = nIt;
asg.acc = acc;
%asg.obj = obj;

prOut;
