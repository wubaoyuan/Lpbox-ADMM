function varargout = oness(varargin)
% Produce numerous one matrix.
%
% Input
%   varargin   -  dimensions
%
% Output
%   varargout  -  numerous ones matrix
%
% History
%   create     -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify     -  Feng Zhou (zhfe99@gmail.com), 09-08-2012

% dimension
for i = 1 : nargout
    varargout{i} = ones(varargin{:});
end