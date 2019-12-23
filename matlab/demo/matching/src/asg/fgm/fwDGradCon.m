function Gr = fwDGradCon(X)
% Compute the gradient of constant function part in Frank-Wolfe algorithm.
%
% Input
%   X       -  correspondence, n1 x n2
%           
% Output    
%   Gr      -  gradient, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-16-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% global variable
global GHHQQG;

Gr = 2 * GHHQQG * X;
