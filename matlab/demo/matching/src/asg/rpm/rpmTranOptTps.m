function [c, d] = rpmTranOptTps(x, y, lamda1, lamda2)
% Compute optimal transformation of TPS.
%
% Input
%   x       -  1st point set, n x dim
%   y       -  2nd point set, m x dim
%   lamda1  -  1st weight
%   lamda2  -  2nd weight
%
% Output
%   c       -  TPS weight, n x 3
%   d       -  affine transformation, 3 x 3
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-05-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 04-23-2012

% homogeneous coordinates
[n, dim] = size(x); 
x = [ones(n, 1), x];

[m, ~] = size(y); 
y = [ones(n, 1), y];

% kernel
K = ctps_gen_K(x, x);

% QR decomposition
[q1, q2, R] = ctps_gen_qr(x);

% solution
[c, d] = ctps_gen_cd_regularized(lamda1, lamda2, q1, q2, R, K, y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q1, q2, R] = ctps_gen_qr(x)
% Genrate QR decomposition of x.
%
% Math
%   x = [q1; q2] * r = q1 * R
%
% Input
%   x   -  point set, n x M
%
% Output
%   q1  -  q1 matrix, n x M
%   q2  -  q2 matrix, n x (n - M)
%   R   -  R matrix, M x M

% dimension
[n, M] = size(x);

[q, r] = qr(x);
q1 = q(:, 1 : M);
q2 = q(:, M + 1 : n);
R = r(1 : M, 1 : M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = ctps_gen_K(x, z)
% Generate TPS kernel matrix.
%
% Input
%   x   -  1st point set, n x M 
%   z   -  2nd point set, m x M
%
% Output
%   K   -  TPS kernel, n x m

% dimension
[n, M] = size(x);
[m, M] = size(z);
dim = M - 1;

% calc. the K matrix.
% 2D: K = r^2 * log r
% 3D: K = -r
K = zeros(n, m);
for it_dim = 1 : dim
    tmp = x(:, it_dim + 1) * ones(1, m) - ones(n, 1) * z(:, it_dim + 1)';
    K = K + tmp .* 2;
end

if dim == 2
    mask = K < 1e-10; % to avoid singularity.
    K = 0.5 * K .* log(K + mask) .* (K > 1e-10);
else
    K = -sqrt(K);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c, d] = ctps_gen_cd_regularized(lamda1, lamda2, q1, q2, R, K, y)
% Calc. regularized TPS c,d.
% (Regularize the affine transformation as well).

% dimension
[n, M] = size(y);
dim = M - 1;

gamma = (q2' * K * q2 + lamda1 * eye(n - M)) \ q2' * y;
c = q2 * gamma;

% add regularization for "d" as well:
% d = inv(R) * q1' * (y-K*q2*gamma);
A = (R' * R + lamda2 * eye(length(R))) \ (R' * q1' * (y - K * q2 * gamma) - R' * R);
d = A + eye(dim + 1, dim + 1);
