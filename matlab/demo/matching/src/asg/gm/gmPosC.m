function X = gmPosC(K, Ct, X0, par)
% Computer a better continuous assingment matrix by refining an old one.
%
% Remark
%   nn = n1 x n2
%
% Input
%   K       -  affinity matrix, [] | nn x nn (sparse)
%   Ct      -  constraint, n1 x n2
%   X0      -  initial assignment, n1 x n2
%   par     -  parameter
%     alg   -  method, 'none' | 'grad' | 'rrwm'
%                'none' : do nothing
%                'grad' : graduate assignment
%                'rrwm' : reweighted random walk matching
%                'quad' : naive quadratic programming
%
% Output
%   X       -  permutation matrix, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 10-15-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-14-2012

% function parameter
alg = par.alg;

% skip
if strcmp(alg, 'none')
    X = X0;
    
% graduate assignment
elseif strcmp(alg, 'grad')
    X = gmPosCGrad(K, X0, par);

% rrwm
elseif strcmp(alg, 'rrwm')
    X = gmPosCRrwm(K, X0, par);

% quad
elseif strcmp(alg, 'quad')
    X = gmPosCQuad(-K, X0, par);

else
    error('unknown method: %s', alg);
end
