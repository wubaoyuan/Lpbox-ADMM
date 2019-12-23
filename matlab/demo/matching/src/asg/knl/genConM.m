function Cs = genConM(Ps)
% Generate confusion matrices for MDA.
%
% Input
%   Ps      -  permutation matrices, k x k x n x n
%
% Output
%   Cs      -  confusion matrices, k x k x n x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-02-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

k = size(Ps, 1);
n = size(Ps, 3);

Cs = zeros(k, k, n, n);

for i = 1 : n
    for j = i + 1 : n
        P = Ps(:, :, i, j);
        
        a = rand(1, k);
        b = float2block(a, 10);
        c = diag(b);
        
        d = rand(k, k);
        d = float2block(d, 3);
        c = c + d;

        Cs(:, :, i, j) = c * P;

        Cs(:, :, j, i) = Cs(:, :, i, j)';
    end
end
