function Bs = matrix2block(n, nB)
% Divide the matrix into almost equal-size blocks.
%
% Input
%   n       -  matrix size
%   nB      -  number of blocks
%
% Output
%   Bs      -  blocks, 1 x nB (cell)
%              each element contains subscript for non-zero elements, 2 x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

Bs = cell(1, nB);

% upper triangular part of matrix
A = ones(n, n);
Au = triu(A);
ind = find(Au);
[I, J] = ind2sub(size(A), ind);
sub = [I'; J']; m = size(sub, 2);

% divide into blocks
vec = 1 : m;
vecBs = vec2block(vec, nB);
for iB = 1 : nB
    vecB = vecBs{iB};
    Bs{iB} = sub(:, vecB); 
end
