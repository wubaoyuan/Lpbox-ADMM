function l = G2L(G)
% Convert class indicator matrix to class label vector.
%
% Example
%   input   -  G = [0 1 0 0; ...      
%                   1 0 0 0; ...  
%                   0 0 1 0]
%   usage   -  l = G2L(G)
%   output  -  l = [2 1 3 0]
%
% Input
%   G       -  class indicator matrix, k x n
%
% Output
%   l       -  class label vector, 1 x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-04-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

% dimension
[k, n] = size(G);

[idxI, idxJ] = find(G);
l = zeros(1, n);
l(idxJ) = idxI;
