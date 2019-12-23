function ind = strloc(strs, list, field)
% Locate the position of the first same string in a string cell array or a struct array.
% with the known field name.
%
% Input
%   strs    -  string, 1 x m (cell)
%   list    -  string cell array or struct array that will be searched, 1 x n (cell) | 1 x n (struct)
%   field   -  field name, can be empty
%              If not empty, list is a struct array.
%
% Output
%   ind     -  position, 1 x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

if ~iscell(strs)
    strs = {strs};
end

% list is a struct array
if ~iscell(list)
    if ~exist('field', 'var')
        error('not enough input');
    end
    list = {list.(field)};
end

% dimension
m = length(strs);
n = length(list);
ind = zeros(1, m);

for i = 1 : m
    for j = 1 : n
        if strcmp(strs{i}, list{j})
            ind(i) = j;
            break;
        end
    end
end
