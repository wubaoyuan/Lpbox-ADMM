function [acc, cost] = evalMDA(Cs, Ps)
% Measure the accuracy of the MDA problem.
%
% Input
%   Cs      -  confusion matrices, k x k x n x n
%   Ps      -  permutation matrices, k x k x n x n
%
% Output
%   acc     -  accuracy
%   cost    -  cost
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-02-2008
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

n = size(Cs, 3);

cost = 0; all = 0;
for i = 1 : n
    for j = i + 1 : n
        C = Cs(:, :, i, j);
        P = Ps(:, :, i, j);

        cost = cost + trace(C * P');
        all = all + sum(C(:));
    end   
end

acc = cost / all;
