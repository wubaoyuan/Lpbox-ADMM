function [idxj, idxs] = divCo(n, m, j)
% Access the count number which has been stored in the specified path.
%
% Input
%   n       -  #total number
%   m       -  #parts
%   j       -  part id
%
% Output
%   idxj    -  index for part j, 1 x nj
%   idxs    -  index set, 1 x m (cell), 1 x ni
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-03-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-16-2012

% divide
ni = ceil(n / m);
idxs = cell(1, m);
for i = 1 : m
    if i < m
        idxs{i} = (i - 1) * ni + 1 : i * ni;
    else
        idxs{i} = (i - 1) * ni + 1 : n;
    end
end

idxj = idxs{j};