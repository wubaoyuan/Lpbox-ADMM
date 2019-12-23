function [G, H, U, V] = gphEg2IncA(Eg, n)
% Obtain incidence matrix for asymmetric edges.
%
% Input
%   Eg      -  graph edge, 2 x (2m)
%   n       -  #nodes
%
% Output
%   G       -  node-edge incidence (starting), n x (2m)
%   H       -  node-edge incidence (ending), n x (2m)
%   U       -  incidence component, n x m
%   V       -  incidence component, n x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-11-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 02-01-2012

% dimension
m = round(size(Eg, 2) / 2);

% incidence matrix
[U, V] = zeross(n, m);
for c = 1 : m
    U(Eg(1, c), c) = 1;
    V(Eg(2, c), c) = 1;
end

% incidence
G = [U, V];
H = [V, U];
