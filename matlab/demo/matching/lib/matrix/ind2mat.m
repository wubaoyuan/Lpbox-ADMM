function A = ind2mat(Ind)
% Convert an index matrix to a sparse matrix.
%
% Input
%   Ind     -  index, 2 x (len + 1)
%
% Output
%   A       -  sparse matrix, m x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 10-30-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-16-2011

% dimension
m = Ind(1, end);
n = Ind(2, end);

% index
idx = sub2ind([m, n], Ind(1, 1 : end - 1), Ind(2, 1 : end - 1));

% store
A = zeros(m, n);
A(idx) = 1;
