function G = L2G(l, k)
% Convert class label vector to class indicator matrix.
%
% Sample
%   input   -  l = [2 1 3 0], k = 3
%   call    -  G = L2G(l, k)
%   output  -  G = [0 1 0 0; ...      
%                   1 0 0 0; ...  
%                   0 0 1 0]
%
% Input
%   l       -  class label vector, 1 x n
%   k       -  #class, use the max(l) if not indicated
%
% Output
%   G       -  class indicator matrix, k x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-04-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

n = length(l);
if ~exist('k', 'var')
    k = max(l);
end

G = zeros(k, n);
idxL = find(l);
idx = sub2ind([k, n], l(idxL), idxL);
G(idx) = 1;

