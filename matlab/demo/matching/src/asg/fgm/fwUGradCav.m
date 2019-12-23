function Gr = fwUGradCav(X)
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
global L KQ;
global HH1 HH2 UU VV UUHH VVHH HUH HVH GKGP;
global IndG1 IndG2 IndG1T IndG2T IndH1 IndH2 IndH1T IndH2T;
global GXG HXH;
global isGXG isHXH;

if isGXG == 0
    GXG = multGXH(IndG1T, X, IndG2);
    isGXG = 1;
end

if isHXH == 0
    HXH = multGXH(IndH1T, X, IndH2);
    isHXH = 1;
end

% gradient: 2 * H1 * ((H1' * X * H2) .* L) * H2' - H1 * (HH1 .* UU) * H1' * X - X * H2 * (HH2 .* VV) * H2';
Gr = 2 * multGXH(IndG1, GXG .* KQ, IndG2T) + GKGP;