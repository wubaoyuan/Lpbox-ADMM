function [idxs, m] = divCo2(n, ni)
% Access the count number which has been stored in the specified path.
%
% Input
%   n       -  #total number
%   ni      -  part size
%
% Output
%   idx     -  index, 1 x m (cell), 1 x ni
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-03-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-30-2012

% divide
m = ceil(n / ni);
idxs = cell(1, m);
for i = 1 : m
    if i < m
        idxs{i} = (i - 1) * ni + 1 : i * ni;
    else
        idxs{i} = (i - 1) * ni + 1 : n;
    end
end
