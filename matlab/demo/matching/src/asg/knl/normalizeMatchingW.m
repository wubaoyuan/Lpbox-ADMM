function [W, D1, D2] = normalizeMatchingW(W, E12, nbIter)
% Normalize edge affinity.
%
% Reference
%   T. Cour, P. Srinivasan, and J. Shi, "Balanced graph matching", In NIPS, 2006
%
% Input
%   W       -  must be tril (use W = tril(W))
%   E12     -  n1 * n2 correspondance matrix.
%              assumes that [I12(:, 1), I12(:, 2)] = find(E12); will be such that: the
%              kth dimensions in W corresponds to the match (I12(k,1),I12(k,2))
%   nbIter  -  #iteration
%
% Output
%
% History
%   create  -  Timothee Cour (timothee.cour@gmail.com), 04-21-2008
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-16-2011

% function option
if nargin < 3
    nbIter = 10;
end
options.normalization = 'iterative';
options.nbIter = nbIter;

E12 = full(E12);
[I12(:, 1), I12(:, 2)] = find(E12);

[n1, n2] = size(E12);
options.n1n2 = [n1, n2];
[W, D1, D2] = mex_normalizeMatchingW(W, I12, options);

W = W / max(max(W));
