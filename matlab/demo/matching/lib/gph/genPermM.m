function Ps = genPermM(k, n)
% Generate Permutation matrices for MDA.
%
% Input
%   k       -  #classes
%   n       -  #subjects
%
% Output
%   Ps      -  permutation matrices, k x k x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-02-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

Ps = zeros(k, k, n);
for i = 1 : n
    Ps(:, :, i) = genPerm(k);
end
