function xs = linspace2(a, b, n, c)
% Decompose a set of n points into m subsets according the specified way.
% The size of subset must sum to n. i.e. sum(ns) == n.
%
% Example
%   input   -  n = 10; m = 3; alg = 'eq';
%   call    -  [ns, idx] = divN(n, m, 'alg', alg)
%   output  -  ns = [3; 3; 4];
%              idx = [2; 1; 2; 3; 3; 3; 1; 2; 3; 1];
%
% Input
%   n       -  number
%   m       -  parts
%   varargin
%     alg   -  algorithm type, 'rand' | 'w' | {'eq'}
%
% Output
%   ns      -  size of subset, m x 1
%   idx     -  index of subset that each number (1 ~ n) has been assigned to, n x 1
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 05-30-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-24-2012

xs = linspace(a, b, n);
ds = abs(xs - c);
[~, p] = min(ds);
xs(p) = c;