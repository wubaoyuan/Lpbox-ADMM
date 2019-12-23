function value = ps(option, name, default, varargin)
% Parse the parameter specified in a struct or in a cell array.
%
% Example1 (when option is a struct):
%   input    -  option.lastname = 'zhou';
%   call     -  value = ps(option, 'lastname', 'noname');
%   output   -  value = 'zhou'
%
% Example2 (when option is a cell array):
%   input    -  option = {'lastname', 'zhou'};
%   call     -  value = ps(option, 'lastname', 'noname');
%   output   -  value = 'zhou'
%
% Input
%   option   -  struct or cell array
%   name     -  field name
%   default  -  default field value
%   values   -  all feasible valids
%
% Output
%   value    -  field value
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

if iscell(option)
    if isempty(option)
        option = [];
    elseif length(option) == 1
        option = option{1};
    else
        option = cell2option(option);
    end
end

if isfield(option, name)
    value = option.(name);
else
    value = default;
end
