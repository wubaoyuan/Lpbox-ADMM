%% Check the answer
% by Jungmin Lee in CVL
function [ nTrueMatch objValue objInValue] = checkAnswer( X, M, T12, Algorithm)

[n1 n2] = size(X);
Tzero = zeros(n1, n2);
Tones = ones(n1, n2);

if sum(sum((X == Tzero) | (X == Tones))) ~= n1*n2;
    X
    error(['Algorithm: ' Algorithm.strName ' returns wrong solution. It is not binary']);
end

from1to2 = sum(X,2);
from2to1 = sum(X,1);

if ~isempty(find(from1to2 > 1)) || ~isempty(find(from2to1 > 1))
    X
    error(['Algorithm: ' Algorithm.strName ' returns wrong solution. It is not 1-to-1']);
end

nTrueMatch = sum(sum(X.*T12));
objValue = X(:)'*M*X(:);
X2 = X.*T12;
objInValue = X2(:)'*M*X2(:);