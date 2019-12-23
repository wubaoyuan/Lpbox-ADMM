function acc = matchAsg(X, asgT)
% Find the match between two assignments and compute the accuracy.
%
% Input
%   X       -  original assignment matrix, n1 x n2
%   asgT    -  ground-truth assignment (can be [])
%
% Output
%   acc     -  accuracy of the matching
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-22-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 02-14-2012

% ground-truth
% XTT = ps(asgT, 'X', []);
% XT = XTT(1:15, 1:15); 
XT = ps(asgT, 'X', []);
if isempty(XT)
    acc = 0;
    return;
end

% non-zero correspondence
idx = find(XT);

% #correct correspondences found
co = length(find(XT(idx) == X(idx)));

% percentage
acc = co / length(idx);
