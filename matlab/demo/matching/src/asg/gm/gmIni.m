function X = gmIni(K, Ct, par)
% Initialize the assingment matrix.
%
% Remark
%   nn = n1 x n2
%
% Input
%   K       -  affinity matrix, nn x nn (sparse)
%   Ct      -  correspondence constraint, n1 x n2
%   par     -  parameter
%     alg   -  method, 'unif' | 'sm' | 'smac'
%                'unif' : uniform value
%                'sm'   : spectral matching (top eigen-vector of K)
%                'smac' : spectral matching with constraint
%
% Output
%   X       -  permutation matrix, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 10-15-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-16-2012

% function parameter
alg = par.alg;

% uniform
if strcmp(alg, 'unif')
    X = gmIniUnif(Ct, par);

% spectral matching    
elseif strcmp(alg, 'sm')
    X = gmIniSm(K, Ct, par);
    
% spectral matching with affine constraint
elseif strcmp(alg, 'smac')
    X = gmIniSmac(K, Ct, par);
    
else
    error('unknown initialization method: %s', alg);
end
