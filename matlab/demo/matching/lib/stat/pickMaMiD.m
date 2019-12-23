function [X, i, ds] = pickMaMiD(Xs)
% Pick the sample set whose minmum pairwise sample distance is maxmimum.
%
% Input
%   Xs      -  sample set, 1 x m (cell), d x n
%
% Output
%   X       -  sample set, d x n
%   i       -  index of the optimum set
%   ds      -  all distances, 1 x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 07-17-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
m = length(Xs);
n = size(Xs{1}, 2);

% distance
ds = zeros(1, m);
for i = 1 : m
    X = Xs{i};

    % pairwise sample distance
    D = conDst(X, X);
    
    % the minimum distance
    D = D + diag(ones(1, n) * inf);
    ds(i) = min(min(D(:)));
end

% the set with maximum distance
[~, i] = max(ds);
X = Xs{i};
