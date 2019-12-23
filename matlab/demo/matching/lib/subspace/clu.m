function G = clu(X, K, k, par)
% An interface of differenct clustering algorithm.
%
% Input
%   X       -  feature matrix, [] | d x n
%   K       -  kernel matrix, [] | n x n
%                if K == [], then a KG has been specified in the global environment
%   k       -  cluster number
%   par     -  parameter
%     alg   -  algorithm name, {'sc'} | 'kmeans'
%
% Output
%   G       -  class indicator matrix, k x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-22-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-25-2012

% function parameter
alg = ps(par, 'alg', 'sc');

if strcmp(alg, 'sc')
    G = cluSc(K, k);
elseif strcmp(alg, 'kmeans')
    [~, G] = kmeanFast(X, k);
else
    error('unknonw alg: %s', alg);
end
