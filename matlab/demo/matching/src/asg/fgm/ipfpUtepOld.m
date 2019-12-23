function X = ipfpUtepOld(K, X0)
% This function tries to maximize the matching score x' K x 
% where x obeys discrete one-to-one matching constraints
% such that x(i) = 1 if nodes(i) is matched to labels(i) and 0 otherwise.
%
% Reference
%   M. Leordeanu and M. Hebert and R. Sukthankar, "An Integer Projected
%   Fixed Point Method for Graph Matching and MAP Inference", in NIPS, 2009
%
% Remark
%   nn = n1 x n2
%
% Input
%   K        -  affinity matrix, nn x nn (sparse)
%   X0       -  initial assignment, n1 x n2
%   par      -  parameter
%     nItMa  -  #maximum iteration steps, {50}
%
% Output
%   X        -  correspondence matrix, n1 x n2
%
% History
%   create  -  Marius Leordeanu (leordeanu@gmail.com), 02-25-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% dimension
[n1, n2] = size(X0);

x0 = X0(:);
a = K * x0;
A = reshape(a, [n1 n2]);

% gradient direction
YY = gmPosDHun(A);
Y = YY - X0;
y = Y(:);

% step size
k = y' * K * y;
if k >= 0
    t = 1;
else
    c = x0' * K * y(:);
    t = min([1, -c / k]);
end
X = X0 + t * Y;
