function X = gmPosCQuad(K, X0, par)
% Quadratic programming for graduate assignment.
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
%   K        -  affinity matrix, nn x nn (sparse)
%   X0       -  initial assignment, n1 x n2
%   par      -  parameter
%     qp     -  toolbox used for quadartic programming, {'mosek'} | 'matlab' | 'cvx' | ...
%                See function qprog for more details
%     b0     -  max(n1, n2) | .5
%     bStep  -  {1.075}
%     bMax   -  {200} | 10
%     tolB   -  {1e-3}
%     tolC   -  {1e-3}
%
% Output
%   X        -  permutation matrix, n1 x n2
%
% History
%   create   -  Timothee Cour (timothee.cour@gmail.com), 02-01-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-22-2011

% dimension
[n1, n2] = size(X0);

% function parameter
qp = ps(par, 'qp', 'matlab');

% A and b
A1 = kron(ones(1, n2), eye(n1));
A2 = kron(eye(n2), ones(1, n1));
A = [A1; A2];
b = ones(n1 + n2, 1);
blx = zeros(n1 * n2, 1);
bux = ones(n1 * n2, 1);

% quadratic programming
x = optQuad(qp, K, [], A, [], b, blx, bux, X0(:));

% vector -> matrix
X = reshape(x, n1, n2);
