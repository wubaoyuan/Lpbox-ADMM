function varargout = cellss(varargin)
% Produce numerous cell matrices.
%
% Example
%   input      -  m = 10; n = 5;
%   call       -  [a, b, c] = cellss(m, n)
%   output     -  a = cell(m, n);
%                 b = cell(m, n);
%                 c = cell(m, n);
%
% Input
%   varargin   -  dimensions
%
% Output
%   varargout  -  cell matrices
%
% History
%   create     -  Feng Zhou (zhfe99@gmail.com), 02-13-2009
%   modify     -  Feng Zhou (zhfe99@gmail.com), 09-08-2012
 
% dimension
for i = 1 : nargout
    varargout{i} = cell(varargin{:});
end
