function As = cellDec(A, s)
% Decompose the matrix into small ones.
%
% Example
%   Input   -  A = [1 3 4 3 1 4; ...   s = [1 5 7];
%                   1 3 3 2 4 5]
%   Call    -  As = cellDec(A, s)
%   Output  -  As = {[1 3 4 3; 1 3 3 2], [1 4; 4 5]}
%
% Input
%   A       -  original matrix, d1 x d2 x ... x dk x n
%   ns      -  segmentation, 1 x (m + 1), d1 x d2 x ... x dk x ni
%
% Output
%   As      -  part set, 1 x m (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-06-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
ndim = length(size(A));
n = size(A, ndim);
ds = cell(1, ndim);
for idim = 1 : ndim
    ds{idim} = size(A, idim);
end
m = length(s) - 1;

As = cell(1, m);
B = reshape(A, [], n);
for i = 1 : m
    As{i} = reshape(B(:, s(i) : s(i + 1) - 1), ds{1 : ndim - 1}, []);
end
