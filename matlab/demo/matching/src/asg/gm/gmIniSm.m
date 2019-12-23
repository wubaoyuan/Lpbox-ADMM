function X = gmIniSm(K, Ct, par)
% Compute the assingment matrix by the algorithm of spectral matching.
%
% Reference
%   M. Leordeanu and M. Hebert, "A Spectral Technique
%   for Correspondence Problems Using Pairwise Constraints", in ICCV, 2005
%
% Math
%   This algorithm is to obtain the optimal x for the following problem
%     max_x   x' * K * x
%     s.t.    x' * x = 1
%
% Remark
%   nn = n1 x n2
%
% Input
%   K       -  affinity matrix, nn x nn (sparse)
%   Ct      -  constraint, n1 x n2
%   par     -  parameter
%     top   -  method of computing the top eigenvector, {'eigs'} | 'pow'
%                'eigs' : the eigs function by MATLAB
%                'pow'  : power iteration
%     nIt   -  iteration number, {30}
%                only used if top == 'pow'
%
% Output
%   X       -  permutation matrix, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-25-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-16-2012

% function parameter
top = ps(par, 'top', 'eigs');
nIt = ps(par, 'nIt', 30);

% dimension
[n1, n2] = size(Ct);
ns = [n1, n2];
nn = ns(1) * ns(2);

% MATLAB
if strcmp(top, 'eigs')
    [x, ~] = eigs(K, 1);
    
% power iteration
elseif strcmp(top, 'pow')
    x = ones(nn, 1);
    x = x / norm(x);
    for iIt = 1 : nIt
        x = K * x;
        x = x / norm(x);
    end

else
    error('unknown method: %s', top);
end

% make sure the eigenvector is positive
% [~, idx] = max(abs(x));
% if x(idx) <= 0
%     x = -x;
% end
x = abs(x);

% vector -> matrix
X = reshape(x, ns);
