function [vals, idxs] = sort2D(A, type)
% Sort the elements of A.
%
% Input
%   A       -  matrix, m x n
%   type    -  sort type, {'ascend'} | 'descend'
%
% Output
%   vals    -  sorted elements, 1 x (m x n)
%   idxs    -  index of element, 2 x (m x n)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-09-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

if ~exist('type', 'var')
    type = 'ascend';
end

[m, n] = size(A); siz = [m, n];
val0s = A(:);
[vals, ind] = sort(val0s, type);

idxs = zeros(2, m * n);
for i = 1 : m * n
    [subi, subj] = ind2sub(siz, ind(i));
    
    idxs(:, i) = [subi; subj];
end
