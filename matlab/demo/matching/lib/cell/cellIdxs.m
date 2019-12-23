function varargout = cellIdxs(idx, varargin)
% Obtain element from a bunch of cell arrays.
%
% Input
%   idx        -  position, 1 x k
%   varargin   -  field name list, 1 x m (cell)
%
% Output
%   varargout  -  field value list, 1 x m (cell), 1 x k (cell)
%
% History
%   create     -  Feng Zhou (zhfe99@gmail.com), 02-01-2009
%   modify     -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

m = nargout;
for i = 1 : m
    varargout{i} = varargin{i}(idx);
end
