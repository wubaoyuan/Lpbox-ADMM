function [G, acc, C, P] = match(G0, varargin)
% Find the best matching of the indicator matrix given the ground-truth label.
%
% Remark
%   The matching can be completed based on the standard indicator matrix
%   or pre-computed confusion matrix.  
%
% Input
%   G0      -  original indicator matrix, k x n
%   varargin
%     GT    -  the ground-truth frame indicator matrix, {[]} | k x n
%     C0    -  confusion matrix before matching, {[]} | k x k
%     P     -  Permuation matrix, {[]} | k x k
%
% Output
%   G       -  new indicator matrix, k x n
%   acc     -  accuracy of the matching
%   C       -  confusion matrix after matching, k x k
%   P       -  Permutation matrix, k x k 
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-22-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
GT = ps(varargin, 'GT', []);
C0 = ps(varargin, 'C0', []);
P = ps(varargin, 'P', []);

% confusion matrix
if ~isempty(GT)
    C0 = GT * G0';
    [P, acc] = asg(C0);

elseif ~isempty(C0);
    [P, acc] = asg(C0);

elseif ~isempty(P)
    C0 = [];
    acc = 0;

else
    error('not enough inputs.');
end
C = C0 * P';

% adjust indicator matrix
k = size(P, 1);
idx = zeros(1, k);
for i = 1 : k
    p = P(i, :);
    idx(i) = find(p);
end
G = G0(idx, :);
