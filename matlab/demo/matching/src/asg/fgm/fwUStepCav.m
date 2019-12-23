function [a, b] = fwUStepCav(X, Y)
% Compute the step size of original function part in Frank-Wolfe algorithm.
%
% Input
%   X       -  correspondence, n1 x n2
%   Y       -  optimal search direction, n1 x n2
%           
% Output    
%   a       -  second-order coefficient
%   b       -  first-order coefficient
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

% auxiliary variables
GYG = multGXH(IndG1T, Y, IndG2);
H1TY = multGXH(IndH1T, Y, []);
YH2 = multGXH([], Y, IndH2);
H1TX = multGXH(IndH1T, X, []);
XH2 = multGXH([], X, IndH2);
HYH = multGXH([], H1TY, IndH2);

% second-order part
a = multTr(KQ .* GYG .^ 2);

% first-order part
tmp1 = multTr(KQ .* GXG .* GYG);
tmp2 = multTr(GKGP .* Y);
b = 2 * tmp1 + tmp2;