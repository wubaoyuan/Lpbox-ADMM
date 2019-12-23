function varargout = matFld(path, varargin)
% Load the mat file from the specified path.
%
% Input
%   path       -  mat path
%   varargin   -  name list, 1 x n (cell)
%
% Output
%   varargout  -  value list, 1 x n (cell)
%
% History
%   create     -  Feng Zhou (zhfe99@gmail.com), 01-06-2009
%   modify     -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

a = load(path);
m = length(varargin);

if nargout ~= m
    error('The number of outputs must be same as inputs.');
end

for i = 1 : m
    varargout{i} = a.(varargin{i});
end
