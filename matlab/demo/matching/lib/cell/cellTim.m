function as = cellTim(varargin)
% Matrix production.
%   as{i, j} = a1s{i, j} * a2s{i, j} * a3s{i, j} ...
%
% Example
%   input     -  a1s = {[1 2 3], [2 4 5]}; a2s = {[.1 .2 .3]', [.1 .2 .3]'};
%   call      -  as = cellTim(a1s, a2s);
%   output    -  as = {[1.4], [2.5]};
%
% Input
%   varargin  -  a set of cell matrices
%
% Output
%   as        -  cell matrix
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 09-20-2010
%   modify    -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

nargin = length(varargin);
as = varargin{1};
m = length(as); 

for i = 2 : nargin
    bs = varargin{i};

    for j = 1 : m
        as{j} = as{j} * bs{j};
    end
end
