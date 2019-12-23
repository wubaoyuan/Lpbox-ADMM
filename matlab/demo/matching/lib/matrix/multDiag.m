function A = multDiag(alg, A0, v)
% Matrix multiply with a vector.
%
% Math
%   alg == 'row'
%     Each row of A0 will be multiple with v, ie, A(:, i) = A0(:, i) .* v.
%   alg == 'col'
%     Each column of A0 will be multiple with v, ie, A(i, :) = A0(i, :) .* v.
%
% Input
%   alg     -  algorithm type, 'row' | 'col'
%   A0      -  original matrix, m x n
%   v       -  vector
%                alg == 'row': v is a 1 x m vector
%                alg == 'col': v is a 1 x n vector
%
% Output
%   A       -  new matrix, m x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-10-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-08-2012

% dimension
[m, n] = size(A0);

if strcmp(alg, 'row')
    V = repmat(v(:), 1, n);
    
elseif strcmp(alg, 'col')
    V = repmat(v(:)', m, 1);
    
else
    error('unknown algorithm: %s', alg);
end

A = A0 .* V;

