function Y = multGXHSlow(IndG, X, IndH)
% Compute the matrix product G * X * H in a much faster way.
%
% Input
%   IndG    -  2 x (lenG + 1)
%   X       -  nG x mH
%   IndH    -  2 x (lenH + 1)
%
% Output
%   Y       -  mG x nH
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 09-20-2012

if isempty(IndG)
    G = eye(size(X, 1));
else
    G = ind2mat(IndG);
end

if isempty(IndH)
    H = eye(size(X, 2));    
else
    H = ind2mat(IndH);
end

Y = G * X * H;
