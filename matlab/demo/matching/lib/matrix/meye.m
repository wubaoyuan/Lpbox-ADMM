function B = meye(A)
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
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-24-2012

% dimension
[n, m] = size(A);

B = zeros(n, m * n);
for i = 1 : n
    idx = (i - 1) * m + 1 : i * m;
    B(i, idx) = A(i, :);
end
