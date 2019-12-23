function [xKs, hstKs] = gmStatK(n1s, n2s, m1s, m2s, ranKs)
% Compute the statistics of the kernel used in graph matching.
%
% Input
%
% Output
% 
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-20-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-25-2013

% dimension
ms = min([m1s; m2s]);
nns = n1s .* n2s;

% rank of K
vKs = ranKs ./ nns;

nBin = 200;
stK = linspace(0, 2, nBin + 1);
hstKs = x2hst(vKs, stK);
xKs = stK(1 : end - 1) + 1 / nBin;

