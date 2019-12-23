function A = fill0(A0, val)
% Set the elements with zero value to the specified value. 
%
% Input
%   A0      -  original matrix, dim1 x dim2 x ...
%   val     -  replaced value, {1}
%
% Output
%   A       -  new matrix, k0 x k
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-04-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

if ~exist('val', 'var')
    val = 1;
end

% index
vis = abs(A0) < eps;

% set to the value
A = A0;
A(vis) = val;
