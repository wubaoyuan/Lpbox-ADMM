function P = genPerm(k)
% Generate Permutation matrix for BMA.
%
% Input
%   k       -  #classes
%
% Output
%   P       -  permutation matrices, k x k
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-02-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 10-09-2011

ind = randperm(k);

P = zeros(k, k);
for c = 1 : k
    P(c, ind(c)) = 1;
end
