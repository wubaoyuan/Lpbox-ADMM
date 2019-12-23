function X = smpGauss(me, Var, n, varargin)
% Sampling from a Gaussian distribution.
%
% Input
%   me       -  mean, d x 1
%   Var      -  variance, d x d
%   n        -  #samples
%   varargin
%     maMiD  -  flag of maximizing the minimum of pairwise distance, 'y' | {'n'}
%     nRep   -  #repetitions (only used if maMiD = 'y'), {100}
%
% Output
%   X        -  sample matrix, d x n
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 01-22-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
isMaMiD = psY(varargin, 'maMiD', 'n');
nRep = ps(varargin, 'nRep', 100);

% dimension
d = size(me, 1);

% repeat sampling until satisfying the constraint
if ~isMaMiD
    nRep = 1;
end

Xs = cell(1, nRep);
for iRep = 1 : nRep
    X0 = randn(d, n);

    % transformation
    [V, D] = eig(Var);
    D2 = sqrt(D);

    Xs{iRep} = V * D2 * X0 + repmat(me, 1, n);
end

if isMaMiD
    X = pickMaMiD(Xs);
else
    X = Xs{1};
end
