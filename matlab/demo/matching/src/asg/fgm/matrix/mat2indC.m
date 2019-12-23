function indCs = mat2indC(A)
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
%   indCs   -  index, 1 x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 10-30-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-20-2012
 
% dimension
[m, n] = size(A);

% index
indCs = zeros(1, n);
for i = 1 : n
    idx = find(A(:, i));
    if length(idx) == 0
        % skip
    elseif length(idx) == 1
        indCs(i) = idx;
    else
        error('incorrect A');
    end
end
