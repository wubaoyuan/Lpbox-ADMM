function varargout = zeross(varargin)
% Produce numerous zero matrix.
%
% Input
%   varargin   -  dimensions
%
% Output
%   varargout  -  numerous zeros matrix
%
% History
%   create     -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify     -  Feng Zhou (zhfe99@gmail.com), 09-08-2012

% dimension
for i = 1 : nargout
    varargout{i} = zeros(varargin{:});
end