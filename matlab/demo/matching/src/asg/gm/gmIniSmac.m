function X = gmIniSmac(K, Ct, par)
% Compute the assingment matrix by the algorithm of spectral matching with constraint.
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
%   Ct      -  constraint, n1 x n2
%                Ct_ij = 1: i and j can be matched
%                Ct_ij = 0: i and j cannot be matched
%   ns      -  #nodes, 1 x 2
%   par     -  parameter
%
% Output
%   X       -  permutation matrix, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-25-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-16-2012

% dimension
[n1, n2] = size(Ct);
ns = [n1, n2];
nn = ns(1) * ns(2);

% initial
X0 = compute_matching_from_W(K, Ct, 3, 'both', 1);

% orthonormalize X
X = computeXorthonormal(X0, Ct);
