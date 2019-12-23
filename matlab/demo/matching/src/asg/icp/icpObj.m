function obj = icpObj(gphs, X, tran)
% Calcuate the objective for ICP.
%
% Input
%   gphs    -  graphs, 1 x 2 (cell)
%   X       -  assignment matrix, n1 x n2
%   tran    -  transformation
%
% Output
%   obj     -  objective
%
% History   
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-26-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 02-16-2012

% element
P1 = gphs{1}.Pt;
P2 = gphs{2}.Pt;

% dimension
[n1, n2] = size(X);

% transform
P2 = tranRun(P2, tran);

% KP
KP = 2 * P1' * P2 - (P1 .* P1)' * ones(2, n2) - ones(n1, 2) * (P2 .* P2);

% obj
obj = multTr(KP .* X);
