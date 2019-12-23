function v = multTr(varargin)
% Trace of multiplication of several matrices. 
%
% v = multTr(A, B)
% v = trace(A^T * B)
%
% Input
%   varargin  -  input matrix, 1 x n (cell)
%
% Output
%   v         -  value
%
% History
%   create    -  Feng Zhou (zhfe99@gmail.com), 03-10-2011
%   modify    -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

A = varargin{1};
for i = 2 : nargin
    A = A .* varargin{i};
end
v = sum(A(:));
