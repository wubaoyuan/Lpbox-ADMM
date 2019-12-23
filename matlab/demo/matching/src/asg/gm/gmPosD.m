function X = gmPosD(K, Ct, X0, par)
% Computer a discrete assingment matrix by rounding from a continuous one.
%
% Input
%   K       -  affinity matrix, [] | nn x nn (sparse)
%   Ct      -  constraints, n1 x n2
%   X0      -  continuous correspondence, n1 x n2
%   par     -  parameter
%     alg   -  method, 'gre' | 'hun' | 'ipfp'
%                'gre'  : greedy algorithm
%                'hun'  : hungraian algorithm
%                'ipfp' : integer fixed point algorithm
%
% Output
%   X       -  discrete correspondence, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 10-07-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-14-2012

% function parameter
alg = par.alg;

if strcmp(alg, 'gre')
    X = gmPosDGre(X0);
    
elseif strcmp(alg, 'hun')
    X = gmPosDHun(X0, Ct);
    
elseif strcmp(alg, 'ipfp')
    X = gmPosDIpfp(K, Ct, X0, par);

else
    error('unknown method: %s', alg);
end
