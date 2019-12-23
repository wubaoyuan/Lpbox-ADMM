function as = cellRep(a, n)
% Replicate matrix into cell array.
%
% Example
%   input     -  a = [1 2 3]; n = 3;
%   call      -  as = cellRep(a, n);
%   output    -  as = {[1 2 3], [1 2 3], [1 2 3]};
%
% Input
%   a         -  matrix
%   n
%
% Output
%   a         -  single cell array, 1 x n (cell)
%   s         -  start position of each segment, 1 x (m + 1)
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 02-16-2009
%   modify    -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

as = cell(1, n);
for i = 1 : n
    as{i} = a;
end
