function B = mtril(A, n)
% Matrix concatenation.
%
% Input
%   A       -  matrix, n x m
%
% Output
%   B       -  matrix, n x mn
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-27-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-19-2013

B1 = tril(A, -n);
B2 = tril(A, -n - 1);
B = B1 - B2;