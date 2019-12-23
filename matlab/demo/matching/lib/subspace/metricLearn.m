function [A, Comp, Lamb, Dire] = metricLearn(X, H, varargin)
% Metric learning with side-information.
%
% Input
%   X0      -  sample matrix, dim x n
%   H0      -  indicator matrix, k x n
%   varargin
%     type  -  {'diag'} | 'full' | 'xing' | 'diag2'
%
% Output
%   A       -  metric parameterization matrix, dim x dim
%   Comp    -  principal components, dim x n
%   Lamb    -  sorted lambda (eigenvalue), dim x 1
%   Dire    -  principal directions, dim x dim
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-18-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
type = ps(varargin, 'type', 'diag');

prompt('metric learning (%s)...\n', type);
if strcmp(type, 'diag')
    A = metricDiag(X, H);

elseif strcmp(type, 'diag2')
    A = metricDiag2(X, H);    

elseif strcmp(type, 'full')
    A = metricFull(X, H);  

elseif strcmp(type, 'xing')
    A = metricXing(X, H);

else
    error('unknown type');
end

[U, D] = svd(A);
D = sqrt(D);
[Lamb, ind] = sort(sum(D, 2), 'descend');
Dire = U(:, ind);
Comp = Dire' * diag(Lamb) * X;

% cost change
in0 = evalClu(X, H, 'type', 'in');
in = evalClu(Comp, H, 'type', 'in');
out0 = evalClu(X, H, 'type', 'out');
out = evalClu(Comp, H, 'type', 'out');
prompt('After learning metric\n');
prompt('  within-class  error %.2f -> %.2f\n', in0, in);
prompt('  between-class error %.2f -> %.2f\n', out0, out);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = metricDiag(X, H)
% Metric learning with diagonal version.
%
% Input
%   X       -  sample matrix, dim x n
%   H       -  indicator matrix, k x n
%
% Output
%   A       -  metric parameterization matrix (diagonal), dim x dim

[dim, n] = size(X);
S = H' * H;
[f, g] = zeross(dim, 1);

for i = 1 : n
    for j = 1 : n
        d = (X(:, i) - X(:, j)) .^ 2;
        if S(i, j)
            f = f + d;
        else
            g = g + d;
        end
    end
end

x = linprog(f, -g', -1, [], [], zeros(dim, 1));
A = diag(sqrt(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = metricFull(X, H)
% Metric learning with full version.
%
% Input
%   X       -  sample matrix, dim x n
%   H       -  indicator matrix, k x n
%
% Output
%   A       -  metric parameterization matrix (full matrix), dim x dim

[dim, n] = size(X);
k = size(H, 1);

% total scatter
Lt = eye(n) * n - ones(n, n);
st = X * Lt * X';

% within-class scatter
sw = zeros(dim, dim);
for c = 1 : k
    vis = H(c, :) == 1;
    Y = X(:, vis);
    m = size(Y, 2);
    
    Lw = eye(m) * m - ones(m, m);
    sw = sw + Y * Lw * Y';
end
sb = st - sw;

fprintf('Optimization...\n');
% cvx_begin
%     variable A(dim, dim) symmetric;
%     minimize(sum(sum((A * WS) .* WS)))
%     subject to
%         A == semidefinite(dim);
%         sum((sum((A * WD) .* WD))) >= 1;
% cvx_end

cvx_begin
    variable A(dim, dim) symmetric;
    minimize(trace(A * sw))
    subject to
        A == semidefinite(dim);
        trace(A * st) >= 1;
cvx_end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = metricXing(X0, H0)
% Metric learning with Eric Xing's version.
%
% Input
%   X0      -  sample matrix, dim x n
%   H0      -  indicator matrix, k x n
%
% Output
%   A       -  metric parameterization matrix (diagonal), dim x dim

X = X0';
S = H0' * H0;
D = 1 - S;
maxiter = 100;

A = opt(X, S, D, maxiter);
