function varargout = cellDim(As, c)
% Obtain dimension of matrices stored in a cell array.
%
% Example1
%   input      -  As = {A11, A12; ...
%                       A21, A22};
%                 A11 = rand(2, 1, 3);
%                 A12 = rand(1, 3, 1);
%                 A21 = rand(3, 2, 4);
%                 A22 = rand(4, 4, 2);
%                 c = 2;
%   call       -  Dim = cellDim(As, c)
%   output     -  Dim = [1 3; ...
%                        2 4];
%
% Example2
%   input      -  As = {A11, A12; ...
%                       A21, A22};
%                 A11 = rand(2, 1, 3);
%                 A12 = rand(1, 3, 1);
%                 A21 = rand(3, 2, 4);
%                 A22 = rand(4, 4, 2);
%   call       -  [Dim1, Dim2, Dim3] = cellDim(As)
%   output     -  Dim1 = [2 1; ...
%                         3 4];
%                 Dim2 = [1 3; ...
%                         2 4];
%                 Dim3 = [3 1; ...
%                         4 2];
%
% Input
%   As         -  set of matrices, n1 x n2 x n3 x ... (cell)
%   c          -  specific dimension (optional)
%
% Output
%   varargout  -  dimension info, n1 x n2 x n3 x ... (when c is provided) | 1 x m (cell) (when c is not provided)
%
% History
%   create     -  Feng Zhou (zhfe99@gmail.com), 12-29-2008
%   modify     -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
ndims = length(size(As));
dims = cell(1, ndims);
for i = 1 : ndims
    dims{i} = size(As, i);
end
dim = numel(As);

% dimension index
k = length(size(As{1}));
if exist('c', 'var')
    idx = c;
else
    idx = 1 : k;
end
if length(idx) ~= nargout
    error('incorrect output number');
end

% output
for j = 1 : length(idx)
    varargout{j} = zeros(dims{:});
    
    for d = 1 : dim
        varargout{j}(d) = size(As{d}, idx(j));
    end
end
