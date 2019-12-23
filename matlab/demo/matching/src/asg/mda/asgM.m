function [Ps, Cs] = asgM(C0s, alg)
% Multi-Dimensional Assignment (MDA).
%
% Input
%   C0s     -  confusion matrix before adjusting, k x k x n x n
%   alg     -  algorithm for finding the matching, {'hug'} | 'bb'
%              'hug': approximated solution by Hungarain algorithm
%              'bb' : exact solution by Branch-and-Bound algorithm (very expensive in time cost)
%
% Output
%   Ps      -  permutation matrix, k x k x n x n
%   Cs      -  confusion matrix after adjusting, k x k x n x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-30-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
[k, k, n, n] = size(C0s);
[Ps, Cs] = zeross(k, k, n, n);

% two subjects
if n == 2
    C0 = C0s(:, :, 1, 2);
    P = asg(C0);

    Ps(:, :, 1, 2) = P;

% multiple (n > 3) subjects
else
    if strcmp(alg, 'hug')
        Ps = asgMHug(C0s);

    elseif strcmp(alg, 'bb')
        Ps = asgMBB(C0s);

    else
        error(['unknown mda algorithm: ' alg]);
    end
end

% adjust the confusion matrix
for i = 1 : n
    Ps(:, :, i, i) = eye(k);

    for j = i + 1 : n
        C0 = C0s(:, :, i, j);
        P = Ps(:, :, i, j);
        Cs(:, :, i, j) = C0 * P';
        Cs(:, :, j, i) = Cs(:, :, i, j)';
    end
end
