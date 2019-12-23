function [a, b] = fwDStepGm(X, Y)
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
global KP KQ;
global IndG1 IndG2 IndG1T IndG2T IndH1 IndH2 IndH1T IndH2T;
global G1s G2s H1s H2s;
global GXG HXH;
global GXGs HXHs;
global isGXG isHXH;

% if isGXG == 0
%     GXG = multGXH(IndG1T, X, IndG2);
%     isGXG = 1;
% end
% 
% if isHXH == 0
%     HXH = multGXH(IndH1T, X, IndH2);
%     isHXH = 1;
% end
% 
% ti = tic;
% GYG = multGXH(IndG1T, Y, IndG2);
% HYH = multGXH(IndH1T, Y, IndH2);
% a = multTr(KQ .* GYG .* HYH);
% b = multTr(KP .* Y) + multTr(KQ, GXG .* HYH + GYG .* HXH);
% t1 = toc(ti);

% ti = tic;
% Xs = sparse(X);
% Ys = sparse(Y);
GYGs = G1s' * Y * G2s;
HYHs = H1s' * Y * H2s;
a = sum(sum(GYGs .* HYHs .* KQ));
b = sum(sum(KP .* Y)) + sum(sum((GXGs .* HYHs + GYGs .* HXHs) .* KQ));
% t2 = toc(ti);
% fprintf('ratio: %f\n', t1 / t2);
% equal('a', a, full(as));
% equal('b', b, full(bs));

