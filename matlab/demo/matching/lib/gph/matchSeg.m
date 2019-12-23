function [seg, acc, C, C0, P] = matchSeg(seg0, segT, varargin)
% Find the best match between two segmentations.
% The segmentation, 'seg0', will be adjusted in its label.
%
% Input
%   seg0    -  original segmentation
%   segT    -  ground-truth segmentation
%     alg   -  measurment of overlapping, {'con'} | 'dis'
%              'con': continuous
%              'dis': discrete
%
% Output
%   seg     -  new segmentation after adjusted the segment's label
%   acc     -  accuracy of the matching
%   C       -  confusion matrix after matching, k x k
%   C0      -  confusion matrix before matching, k x k
%   P       -  permutation matrix, k x k
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-22-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 02-28-2012

% function option
alg = ps(varargin, 'alg', 'con');

% confusion matrix
G0 = seg0.G;

% confusion matrix (continuous)
if strcmp(alg, 'con')
    C0 = genConMatrix(segT, seg0);

% confusion matrix (discrete)
elseif strcmp(alg, 'dis')
    C0 = genConMatrixD(segT, seg0);

else
    error('unknown alg');
end

% optimal assignment via Hungrian algorithm
P = asgHun(C0, 'opt', 'max');
C = C0 * P';
acc = sum(diag(C)) / sum(C(:));

% adjust indicator matrix
k = size(P, 1);
idx = zeros(1, k);
for i = 1 : k
    p = P(i, :);
    idx(i) = find(p);
end
G = G0(idx, :);

% store
seg = seg0;
seg.G = G;
seg.acc = acc;
