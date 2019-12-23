function Ind = mat2ind(A)
% Convert a sparse binary matrix to an index matrix.
%
% Example
%   input   -  A = [1 0 0; 
%                   0 1 0]
%   call    -  Ind = mat2ind(A)
%   output  -  Ind = [1 2 2; 
%                     1 2 3]
%
% Input
%   A       -  sparse matrix, m x n
%
% Output
%   Ind     -  index, 2 x (len + 1)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 10-30-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 02-19-2012
 
% dimension
[m, n] = size(A);

% index
[is, js] = find(A);
len = length(is);

% store
Ind = [is, js]';
Ind = [Ind, [m; n]];
