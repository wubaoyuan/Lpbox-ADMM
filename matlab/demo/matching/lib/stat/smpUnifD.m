function X = smpUnifD(mi, ma, n, varargin)
% Sampling from a uniform distribution.
% Notice that the sample value is discrete, i.e.,
%   each sample is an integer between minimum and maximum interger.
%
% Input
%   mi      -  minimum of each dimension, d x 1
%   ma      -  maximum of each dimension, d x 1
%   n       -  #samples
%
% Output
%   X       -  sample, d x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 07-19-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

X = genUnif(mi, ma, n);
X = round(X);
