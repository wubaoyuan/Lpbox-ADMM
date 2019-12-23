function str = cellStr(varargin)
% Representing a cell array as a string.
%
% Input
%   varargin  -  set of cell array, 1 x m (cell)
%
% Output
%   str       -  string
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 11-18-2010
%   modify    -  Feng Zhou (zhfe99@gmail.com), 06-15-2012

% dimension
m = nargin;

strs = cell(1, m);
for i = 1 : m
    a = varargin{i};
    if isa(a, 'numeric') || isa(a, 'logical')
        strs{i} = mat2str(a);
    elseif isa(a, 'char')
        strs{i} = a;
    elseif isa(a, 'struct')
        strs{i} = struct2str(a);
    else
        error('unknown type')
    end
    
    if i > 1
        strs{i} = ['_' strs{i}];
    end
end
str = [strs{:}];

%%%%%%%%%%%%%%%%%%%%%%%%%
function str = mat2str(a)

n = numel(a);

strs = cell(1, n);
for i = 1 : n
    strs{i} = num2str(a(i));
    if i > 1
        strs{i} = ['_' strs{i}];
    end
end
str = [strs{:}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = struct2str(a)

% all fields
nms = fieldnames(a);
n = length(nms);

% each field
strs = cell(1, n);
for i = 1 : n
    nm = nms{i};
    val = a.(nm);
    
    if isa(val, 'numeric')
        strs{i} = [nm '-' mat2str(val)];
    elseif isa(val, 'char')
        strs{i} = [nm '-' val];
    elseif isa(val, 'struct')
        strs{i} = [nm '-' struct2str(val)];
    else
        error('unknown type')
    end
    if i > 1
        strs{i} = ['_' strs{i}];
    end
end
str = [strs{:}];