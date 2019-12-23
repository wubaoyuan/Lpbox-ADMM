function X = ipfpUStep(H1, H2, X0)
% This function tries to maximize the matching score x' K x.
%
% Reference
%   M. Leordeanu and M. Hebert and R. Sukthankar, "An Integer Projected
%   Fixed Point Method for Graph Matching and MAP Inference", in NIPS, 2009
%
% Remark
%   The code is equivalent to the original IPFP algorithm.
%   Instead of using the large K, the function requires the smaller, KP and KQ.
%
% Input
%   H1      -  augment node-edge adjacency matrix, n1 x (m1 + n1)
%   H2      -  augment node-edge adjacency matrix, n2 x (m2 + n2)
%   X0      -  initial assignment, n1 x n2
%
% Output
%   X       -  correspondence matrix, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-20-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

global L;

% dimension
[n1, n2] = size(X0);

% A(:) = K * x0
HXH = H1' * X0 * H2;
A = H1 * (HXH .* L) * H2';

% gradient direction
YY = gmPosDHun(A);
Y = YY - X0;
y = Y(:);

HYH = H1' * Y * H2;

% step size (y' * K * y)
k = multTr(L, HYH .^ 2);

if k >= 0
    t = 1;
else
    % x0' * K * y(:);
    c = multTr(L, HXH, HYH);
    t = min([1, -c / k]);
end
X = X0 + t * Y;
