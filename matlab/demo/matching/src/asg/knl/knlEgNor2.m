function K = knlEgNor2(K0, ns)
% Normalize edge affinity.
%
% Reference
%   T. Cour, P. Srinivasan, and J. Shi, "Balanced graph matching", In NIPS, 2006
%
% Input
%   KQ0      -  original edge affinity, m1 x m2
%   par      -  parameter
%     nor    -  flag of whether to normalize edge affinity, 'y' | {'n'}
%     nItMa  -  #maximum iteration, {10}
%     th     -  threshold, {1e-7}
%
% Output
%   KQ       -  new edge affinity, m1 x m2
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 08-09-2011
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-14-2011

E12 = ones(ns);

% for memory efficiency
K = tril(sparse(K0));

% kronecker normalization
K = normalizeMatchingW(K, E12);

% symmetric
K = tril(K, -1)' + K;
