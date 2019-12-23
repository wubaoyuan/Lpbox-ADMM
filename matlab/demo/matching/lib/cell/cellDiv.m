function as = cellDiv(a, s)
% Divide a cell array into a set of cell arrays.
%
% Example
%   input   -  a = {'feng', 28, 'age', 33};
%              s = [1 3 4 5];
%   call    -  as = cellDiv(a, s);
%   output  -  as{1} = {'feng', 28}; as{2} = {'age'}; as{3} = 33;
%
% Input
%   a       -  single cell array, 1 x n (cell)
%   s       -  start position of each segment, 1 x (m + 1)
%
% Input
%   as      -  set of cell array, 1 x m (cell)
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 09-25-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

m = length(s) - 1;
as = cell(1, m);
for i = 1 : m
    as{i} = {a{s(i) : s(i + 1) - 1}};
end
