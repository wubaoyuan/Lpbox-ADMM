function C = ominus(A, B)
% Title minus.
%
% Example
%    C_{6 x 2} = A_{6 x 2} \ominus B_{2 x 2}
%              = A_{6 x 2} - (1_{3 x 1} \otimes B_{2 x 2})
%
% Input
%   A       -  matrix A, m1 x n
%   B       -  matrix B, m2 x n
%
% Output
%   C       -  matrix C, m1 x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 09-20-2010
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
m1 = size(A, 1);
m2 = size(B, 1);

t = round(m1 / m2);
G = ones(t, 1);

C = A - kron(G, B);
