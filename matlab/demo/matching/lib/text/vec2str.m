function [str, strs] = vec2str(vecs, form, varargin)
% Print a vector with the specified type.
%
% Input
%   vecs     -  scalar vectors, 1 x m (cell), 1 x n
%   form     -  printf format, [] | '%d' | '%f'
%   varargin
%     delim  -  seperator for concatenate strings, {''}
%
% Output
%   str      -  string by concatnating together
%   strs     -  string for each value, 1 x n (cell)
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 01-21-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
delim = ps(varargin, 'delim', '');
if strcmp(delim, '\n')
    delim = newLine;
end

if ~iscell(vecs)
    vecs = {vecs};
end
m = length(vecs);
n = length(vecs{1});

[strs, str1s] = cellss(1, n);
for i = 1 : n
    % string without delim
    vars = cell(1, m);
    for j = 1 : m
        vars{j} = vecs{j}(i);
    end
    strs{i} = sprintf(form, vars{:});
    
    if strcmp(strs{i}, '0.0')
        strs{i} = '0';
    elseif length(strs{i}) >= 2 && strcmp(strs{i}(1 : 2), '0.')
        strs{i} = strs{i}(2 : end);
    end
    
    % string with delim
    str1s{i} = strs{i};
    if ~isempty(delim) && i < n
        str1s{i} = [str1s{i}, delim];
    end
end
str = cat(2, str1s{:});
