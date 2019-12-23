function varargout = setIdx(idx, varargin)
% Fetch out the field of each cell element.
%
% Input
%   idx        -  index
%   varargin   -  set of cell array, 1 x m (cell)
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
    if isa(varargin{i}, 'cell')
        varargout{i} = varargin{i}{idx};
    elseif isa(varargin{i}, 'numeric')
        varargout{i} = varargin{i}(idx);
    elseif isa(varargin{i}, 'struct')
        varargout{i} = varargin{i}(idx);
    end
end
