function M = mdiag(M, v)
% Update the diagonal of a matrix.
%
% Input
%   M       -  matrix, n x n
%   v       -  scalar or vector, 1 x 1 | 1 x n
%
% Output
%   M       -  matrix
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-29-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 01-29-2012

% dimension
n = size(M, 1);

% index
idx = sub2ind([n n], 1 : n, 1 : n);

% row vector
if size(v, 1) > 1
    v = v';
end

% set the new value
M(idx) = v;
