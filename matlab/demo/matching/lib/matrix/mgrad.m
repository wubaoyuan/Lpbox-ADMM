function F = mgrad(n)
% Obtain the differential operator used for computing derivative.
%
% Remark
%   The output F is a difference matrix such that:
%     X * F == gradient(X);
%
% Input
%   n       -  #sample
%
% Output
%   F       -  1st order differential operator, n x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 12-18-2008
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

F = zeros(n, n);
idx = n + 1 : n + 1 : n * n - n - 1;
F(idx) = -.5;
F(idx + 2) = .5;

F(1, 1) = -1;
F(2, 1) = 1;
F(end - 1, end) = -1;
F(end, end) = 1;

