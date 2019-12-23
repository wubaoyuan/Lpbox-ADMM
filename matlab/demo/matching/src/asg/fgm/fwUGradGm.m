function Gr = fwUGradGm(X)
% Compute the gradient of original function part in Frank-Wolfe algorithm.
%
% Input
%   X       -  correspondence, n1 x n2
%           
% Output    
%   Gr      -  gradient, n1 x n2
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 02-16-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-06-2013

% global variable
global L;
global IndG1 IndG2 IndG1T IndG2T IndH1 IndH2 IndH1T IndH2T;
global HXH;
global isHXH;

if isHXH == 0
    HXH = multGXH(IndH1T, X, IndH2);
    isHXH = 1;
end

% compute gradient: Gr(:) = K * x0
Gr = multGXH(IndH1, HXH .* L, IndH2T);
