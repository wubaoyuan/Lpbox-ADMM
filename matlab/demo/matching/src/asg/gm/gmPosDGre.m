function X = gmPosDGre(X0)
% Post-processing the continuous assignment to obtain a discrete 
% one by a greedy algorithm.
%
% Reference
%   M. Leordeanu and M. Hebert, "A Spectral Technique
%   for Correspondence Problems Using Pairwise Constraints", in ICCV, 2005
%
% Input
%   X0      -  continuous correspondence, n1 x n2
%
% Output
%   X       -  discrete correspondence, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-25-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-22-2011

% dimension
[n1, n2] = size(X0);
E12 = ones(n1, n2);
X0(E12 == 0) = -1;

% greedy discretization
X = zeros(n1, n2);
[val, idx] = max(X0(:));
n = 0;
while val > 0
    X(idx) = 1;
    n = n + 1;

    % index
    [i, j] = ind2sub([n1 n2], idx);

    % remove
    X0(:, j) = 0;
    X0(i, :) = 0;

    % next position
    [val, idx] = max(X0(:));
end
