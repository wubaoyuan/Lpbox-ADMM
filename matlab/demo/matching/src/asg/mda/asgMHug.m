function Ps = asgMHug(Cs)
% Approximate the MDA solution by using Hungarain algorithm.
%
% Input
%   Cs      -  confusion matrix, k x k x n x n
%
% Output
%   Ps      -  permutation matrix, k x k x n x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-03-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
n = size(Cs, 3);

cost = -inf;
for i = 1 : n
    P2s = asgOne(Cs, i);
    [acc2, cost2] = evalMda(Cs, P2s);

    if cost2 > cost
        cost = cost2;
        Ps = P2s;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ps = asgOne(Cs, r)
% Do assignment with the specific subject.
%
% Input
%   Cs  -  confusion matrices, k x k x n x n
%   r   -  index of reference point
%
% Output
%   Ps  -  permutation matrices, k x k x n x n

k = size(Cs, 1);
n = size(Cs, 3);

Ps = zeros(k, k, n, n);

for j = 1 : n
    C = Cs(:, :, r, j);
    Ps(:, :, r, j)  = asg(C);
    Ps(:, :, j, r) = Ps(:, :, r, j)';
end

for i = 1 : n
    if i == r, continue; end
    for j = i + 1 : n
        if j == r, continue; end
        Ps(:, :, i, j) = Ps(:, :, i, r) * Ps(:, :, r, j);
        Ps(:, :, j, i) = Ps(:, :, i, j)';
    end
end
