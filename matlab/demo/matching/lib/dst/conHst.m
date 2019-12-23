function hst = conHst(X, varargin)
% Obtain a histogram for the given sample space.
%
% Input
%   X       -  sample matrix, dim x n
%   varargin
%     alg   -  algorithm to divide the sample space, {'eq'} | 'kmeans'
%              'eq':     equally divided
%              'kmeans': non-equal divided ('kmeans')
%              'gmm':    non-equal divided ('GMM')
%     val   -  value dependent on the algorithm, {[]}
%
% Output
%   hst    -  histogram
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-11-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% function option
alg = ps(varargin, 'alg', 'eq');
val = ps(varargin, 'val', []);
hst.alg = alg;
hst.val = val;

if strcmp(alg, 'eq')
    hst.mi = min(X, [], 2);
    hst.ma = max(X, [], 2);

elseif strcmp(alg, 'kmeans')
else
    error('unknown method');
end
