function [a, b] = fwUStepGm(X, Y)
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
global L;
global IndG1 IndG2 IndG1T IndG2T IndH1 IndH2 IndH1T IndH2T;
global GXG HXH;
global isGXG isHXH;

if isHXH == 0
    HXH = multGXH(IndH1T, X, IndH2);
    isHXH = 1;
end

HYH = multGXH(IndH1T, Y, IndH2);
a = multTr(L, HYH, HYH);
b = 2 * multTr(L, HXH, HYH);
