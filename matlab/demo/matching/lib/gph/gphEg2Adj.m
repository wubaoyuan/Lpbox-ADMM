function [A, idx] = gphEg2Adj(Eg, n)
% Obtain the adjacency matrix from edge index matrix.
%
% Input
%   Eg      -  graph edge, 2 x 2m
%   n       -  #nodes
%
% Output
%   A       -  node-node adjacency, n x n
%   idx     -  index of each edge, 1 x 2m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 02-26-2013

% dimension
m = round(size(Eg, 2) / 2);

A = zeros(n, n);
idx = sub2ind([n n], Eg(1, :), Eg(2, :));
A(idx) = 1;