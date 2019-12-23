function varargout = mat2inds(varargin)
% Convert a sparse binary matrix to an index matrix.
%
% Input
%   varargin   -  sparse matrices, 1 x k (cell)
%
% Output
%   varargout  -  index matrices, 1 x k (cell)
%
% History
%   create     -  Feng Zhou (zhfe99@gmail.com), 02-15-2012
%   modify     -  Feng Zhou (zhfe99@gmail.com), 02-15-2012

for c = 1 : nargin
    varargout{c} = mat2ind(varargin{c});
end
