function [nP, nZ, nN, ma, mi, gap] = eign(A)
% Calculate the number of positive, zero and negative eigenvalues of A.
%
% Input
%   A       -  input matrix, n x n
%
% Output
%   nP      -  #positive eigenvalues
%   nZ      -  #zero eigenvalues
%   nN      -  #negative eigenvalues
%   ma
%   mi
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 10-18-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-20-2011

% dimension
n = size(A, 1);

% eigen-decomposition
if issparse(A)
    A = full(A);
end
[~, D] = eig(A);    
ds = sort(diag(D), 'descend');

% account
nP = 0;
nZ = 0;
nN = 0;
for i = 1 : n
    if ds(i) > eps
        nP = nP + 1;
    elseif ds(i) < -eps
        nN = nN + 1;
    else
        nZ = nZ + 1;
    end
end

ma = ds(1);
mi = ds(end);

gap = ds(1) - ds(2);
