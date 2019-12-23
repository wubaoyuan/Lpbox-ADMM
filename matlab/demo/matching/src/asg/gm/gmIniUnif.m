function X = gmIniUnif(Ct, par)
% Compute the assingment matrix by uniformly assigning the value of X.
%
% Input
%   Ct      -  constraints, n1 x n2
%   par     -  parameter
%     nor   -  algorithm for normalization, {'none'} | 'unit' | 'doub'
%                'none' : no normalization on X
%                'unit' : unit normalization on vec(X)
%                'doub' : X has to be a doubly stochastic matrix
%
% Output
%   X       -  continuous correspondence matrix, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-25-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-14-2012

% function parameter
nor = ps(par, 'nor', 'none');

% dimension
[n1, n2] = size(Ct);
ns = [n1, n2];

% init
X = ones(ns) + eps;

% constraints
X(Ct == 0) = 0;

% post-processing
if strcmp(nor, 'none')
    
elseif strcmp(nor, 'unit')
    X = X ./ norm(X(:));
    
elseif strcmp(nor, 'doub')
    X = bistocNormalize_slack(X, 1e-7);
    
else
    error('unknown algorithm: %s', non);
end
