function [G, Y] = cluSc(K, k)
% Spectral clustering.
%
% Reference
%   A. Y. Ng, M. I. Jordan and Y. Weiss, "On spectral clustering: analysis and an algorithm," in NIPS, 2001
%
% Input
%   K       -  kernel matrix, [] | n x n
%              if K == [], then a KG has been specified in the global environment
%   k       -  cluster number
%
% Output
%   G       -  indicator matrix, k x n
%   Y       -  sample matrix after embedding, k x n
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-22-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 07-25-2012

% global kernel
global KG;
isKG = isempty(K);

% Formula:
% De = diag(sum(K, 2));
% De2 = sqrt(De);
% L = De2 \ K / De2;
if isKG
    de2 = sqrt(sum(KG, 2));
    L = KG ./ (de2 * de2');
else
    de2 = sqrt(sum(K, 2));
    L = K ./ (de2 * de2');
end

X = eigk(L, k);
Y = X';

% normalize rows of X
tmp = sqrt(sum(X .^ 2, 2));
for i = 1 : length(tmp)
    if abs(tmp(i)) < eps
        tmp(i) = 1;
    end
end
X = diag(tmp) \ X;

% k-means
[~, labs] = kmeanFast(X, k);
G = L2G(labs, k);
