function [idx1, idx2] = asgX2Idx(X)
% Convert correspondence matrix to index.
%
% Input
%   X       -  correspondence matrix, n1 x n2
%
% Output
%   idx1    -  index 1, 1 x n1
%   idx2    -  index 2, 1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-12-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 01-12-2012

% dimension
[n1, n2] = size(X);

idx1 = zeros(1, n1);
for i = 1 : n1
    js = find(X(i, :));
    
    if ~isempty(js)
        idx1(i) = js(1);
    end
end

idx2 = zeros(1, n2);
for j = 1 : n2
    is = find(X(:, j));
    
    if ~isempty(is)
        idx2(j) = is(1);
    end
end
