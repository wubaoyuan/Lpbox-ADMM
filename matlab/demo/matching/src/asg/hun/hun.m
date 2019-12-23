function asg = hun(K, asgT)
% Hungarian algorithm for linear asignment problem.
%
% Math
%   This algorithm is to obtain the optimal X for the following problem
%     max_X   trace(K * X')
%     s.t.    X is a permutation matrix
%
% Input
%   K       -  confusion matrix, n1 x n2
%   asgT    -  ground-truth assignment (can be [])
%
% Output
%   asg     -  asignment
%     alg   -  'hun'
%     X     -  permutation matrix, n1 x n2
%     acc   -  accuracy
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-25-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
[n1, n2] = size(K);
prIn('hun', 'n1 %d, n2 %d', n1, n2);

% convert to minimum-cost assignment problem
K = max(K(:)) - K;
X = hungarian(K);

% compare with ground-truth
acc = matchAsg(X, asgT);

% store
asg.alg = 'hun';
asg.X = X;
asg.acc = acc;

prOut;
