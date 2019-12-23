function X = newAsgRand(n1, n2, n)
% Generate random assignment.
%
% Input
%   n1      -  #points in 1st graph
%   n2      -  #points in 2nd graph
%   n       -  #pairs (can be [])
%
% Output
%   X       -  correspondence matrix, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 09-12-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 02-02-2012

idx1 = randperm(n1);
idx2 = randperm(n2);

if isempty(n)
    n = min(n1, n2);
end

X = zeros(n1, n2);
for i = 1 : n
    X(idx1(i), idx2(i)) = 1;
end
