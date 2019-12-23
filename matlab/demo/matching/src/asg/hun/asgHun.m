function X = asgHun(C, varargin)
% Solve linear assignment by the Hungrian algorithm.
%
% Math
%   X = argmax or argmin trace(C' * X), subject to X is a permutation matrix
%
% Input
%   C       -  weight, n1 x n2
%   varargin
%     opt   -  optimization operator, 'min' | {'max'}
%
% Output
%   X       -  discrete correspondence, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-25-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 02-28-2012

% function option
opt = ps(varargin, 'opt', 'max');

if strcmp(opt, 'max')
    C = max(C(:)) - C;
elseif strcmp(opt, 'min')
    % do noting
else
    error('unknown operator: %s', opt);
end

% c++ implementation
ind2 = assignmentoptimal(C);

% dimension
[n1, n2] = size(C);

% index -> matrix
if n1 <= n2
    idx = sub2ind([n1 n2], 1 : n1, ind2');
    
else
    ind1 = find(ind2);
    ind2 = ind2(ind1);
    idx = sub2ind([n1 n2], ind1', ind2');
end
X = zeros(n1, n2);
X(idx) = 1;

% make sure it is right
%X2 = hungarian(X0);
%equal('X', X, X2);
