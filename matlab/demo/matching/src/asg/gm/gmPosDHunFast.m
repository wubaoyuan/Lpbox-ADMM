function X = gmPosDHunFast(X0, Ct, varargin)
% Post-processing the continuous correspondence matrix
% to obtain a discrete solution by the Hungrian algorithm.
%
% Input
%   X0      -  continuous correspondence, n1 x n2
%   Ct      -  constraint matrix, n1 x n2 | []
%                Con_ij = 1: i and j can be matched
%                Con_ij = 0: i and j cannot be matched
%   varargin
%     opt   -  optimization operator, 'min' | {'max'}
%
% Output
%   X       -  discrete correspondence, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-25-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-26-2012

% function option
opt = ps(varargin, 'opt', 'max');

% dimension
[n1, n2] = size(X0);

idx1 = find(Ct == 1);
X1 = reshape(X0(idx1), [], n1);

% max
if strcmp(opt, 'max')
    [~, idx2] = max(X1);
elseif strcmp(opt, 'min')    
    [~, idx2] = min(X1);    
else
    error('unknown operator: %s', opt);    
end

% index -> matrix
idx3 = sub2ind(size(X1), idx2, 1 : n1);
idx = idx1(idx3);

X = zeros(n1, n2);
X(idx) = 1;
