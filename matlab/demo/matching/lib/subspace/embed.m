function Y = embed(X, parPca)
% Embed samples according to the specific algorithm.
%
% Input
%   alg     -  algorithm type, 'pca' | 'lda' | 'kpca' | 'mds' | 'lle' | 'sc'
%   X       -  original sample matrix, dim0 x n
%   dim     -  dimensions after projection
%
% Output
%   Y       -  new sample matrix, dim x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-25-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011


