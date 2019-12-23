function psCheck(option, varargin)
% Check the option only contains the specific fields.
% If the field does not exist in option, an error will be throwed.
%
% Example:
%   assume    -  option = {'name', 'feng', 'sex', 'male'}
%   call      -  psCheck(option, 'name', 'gender', 'age')
%   then      -  an error will be throwed because the 'sex' field does not exist
%
% Input
%   option    -  name list, 1 x n (cell)
%   varargin  -  name list, 1 x m (cell)
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 03-20-2009
%   modify    -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

m = length(varargin);
for i = 1 : m
    name = varargin{i};
    a.(name) = 'tmp';
end

n = length(option);
for i = 1 : 2 : n
    name = option{i};
    
    if ~isfield(a, name)
        error(['unknown field name: ' name]);
    end
end
