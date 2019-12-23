function varargout = stFld(a, varargin)
% Fetch out the field value according to the specific field name list.
%
% Example
%   input      -  a.name = 'feng';
%                 a.age = 25;
%                 a.gender = 'male'.
%   call       -  [name, age, gender] = stFld(a, 'name', 'age', 'gender')
%   output     -  name = 'feng';
%                 age = 25;
%                 gender = 'male'
%
% Input
%   a          -  struct
%   varargin   -  field name list, 1 x m (cell)
%
% Output
%   varargout  -  field value list, 1 x m (cell)
%
% History
%   create     -  Feng Zhou (zhfe99@gmail.com), 02-01-2009
%   modify     -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

m = length(varargin);

if nargout ~= m
    error('The number of outputs must be same as the names given in the parameters.');
end

for i = 1 : m
    varargout{i} = a.(varargin{i});
end
