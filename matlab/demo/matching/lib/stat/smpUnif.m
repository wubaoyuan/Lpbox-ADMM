function X = smpUnif(mi, ma, n, varargin)
% Sampling from a uniform distribution.
%
% Input
%   mi       -  minimum value of each dimension, d x 1
%   ma       -  maximum value of each dimension, d x 1
%   n        -  #samples
%   varargin
%     maMiD  -  flag of maximizing the minimum of pairwise distance, 'y' | {'n'}
%     nRep   -  #repetitions (only used if maMiD = 'y'), {100}
%
% Output
%   X        -  sample matrix, d x n
%
% History
%   create   -  Feng Zhou (zhfe99@gmail.com), 07-17-2009
%   modify   -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
isMaMiD = psY(varargin, 'maMiD', 'n');
nRep = ps(varargin, 'nRep', 100);

% dimension
d = size(mi, 1);

% repeat sampling until satisfying the constraint
if ~isMaMiD
    nRep = 1;
end

Xs = cell(1, nRep);
for iRep = 1 : nRep
    Xs{iRep} = rand(d, n) .* repmat(ma - mi, 1, n) + repmat(mi, 1, n);
end

if isMaMiD
    X = pickMaMiD(Xs);
else
    X = Xs{1};
end
