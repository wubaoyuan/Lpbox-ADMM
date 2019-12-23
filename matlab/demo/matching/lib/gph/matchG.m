function [G, acc, P, C, C0] = matchG(G0, GT)
% Find the best matching of the indicator matrix given the ground-truth label.
%
% Input
%   G0      -  original indicator matrix, k x n
%   GT      -  ground-truth frame indicator matrix, k x n
%
% Output
%   G       -  new indicator matrix, k x n
%   acc     -  accuracy of the matching
%   P       -  permutation matrix, k x k 
%   C       -  confusion matrix after matching, k x k
%   C0      -  confusion matrix before matching, k x k
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-22-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% confusion matrix
C0 = GT * G0';

% optimal alignment via Hungrian
asg = hun(C0, 'opt', 'max');
[P, acc] = stFld(asg, 'P', 'acc');

C = C0 * P';

% adjust indicator matrix
k = size(P, 1);
idx = zeros(1, k);
for i = 1 : k
    p = P(i, :);
    idx(i) = find(p);
end
G = G0(idx, :);
