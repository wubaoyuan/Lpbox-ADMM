function flag = cluEmp(G)
% Test whether any cluster is empty or not.
%
% Input
%   G       -  indicator matrix, k x n
%
% Output
%   flag    -  empty flag
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

flag = ~isempty(find(sum(G, 2) == 0, 1));
