function [a, s] = cellCat(varargin)
% Concatenate a set of cell array into a long cell array.
%
% Example1
%   input     -  a1 = {'feng', 28}; a2 = {'age'}; a3 = 33;
%   call      -  [a, s] = cellCat(a1, a2, a3);
%   output    -  a = {'feng', 28, 'age', 33};
%                s = [1 3 4 5];
%
% Example2
%   input     -  as = {{'feng', 28}, {'age'}, 33};
%   call      -  [a, s] = cellCat(as);
%   output    -  a = {'feng', 28, 'age', 33};
%                s = [1 3 4 5];
%
% Input
%   varargin  -  set of cell array, 1 x m (cell)
%
% Output
%   a         -  cell array, 1 x n (cell)
%   s         -  start position of each segment, 1 x (m + 1)
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 02-16-2009
%   modify    -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
m = length(varargin);
maxN = 10000;

if m == 1
    as = varargin{1};
    m = length(as);
else
    as = varargin;
end

a = cell(1, maxN); n = 0;
s = ones(1, m + 1);
for i = 1 : m
    a0 = as{i};
    if ~iscell(a0)
        a0 = {a0};
    end

    mi = length(a0);
    for j = 1 : mi
        n = n + 1;
        a{n} = a0{j};
    end
    s(i + 1) = s(i) + mi;
end
a(n + 1 : end) = [];
