function [a, b] = fwDStepCon(X, Y)
% Compute the step size of constant function part in Frank-Wolfe algorithm.
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
global G1 G2 H1 H2 HH1 HH2;
global G1s G2s H1s H2s HH1s HH2s;
global IndG1 IndG2 IndG1T IndG2T IndH1 IndH2 IndH1T IndH2T IndHH1 IndHH2;
global indG1 indG2 indH1 indH2;
global QQ1 QQ2;

X = full(X);
Y = full(Y);
YY = Y * Y';
XY = X * Y';

ti = tic;

% trace(Q1' * Q1 * (G1' * Y * Y' * G1) .* (H1' * H1))
tmp1 = multGXHSQTr(indG1', YY, indG1, IndHH1, QQ1);

% trace(QQ2' * (G2' * Y * Y' * G2) .* (H2' * H2))
tmp2 = multGXHSQTr(indG2', YY, indG2, IndHH2, QQ2);

a = tmp1 + tmp2;

% trace(KQ' * (G1' * X * Y' * G1) .* (H1' * H1))
tmp1 = multGXHSQTr(indG1', XY, indG1, IndHH1, QQ1);

% trace(KQ' * (G2' * X * Y' * G2) .* (H2' * H2))
tmp2 = multGXHSQTr(indG2', XY, indG2, IndHH2, QQ2);

b = 2 * (tmp1 + tmp2);
t1 = toc(ti);

% ti = tic;
% Xs = sparse(X);
% Ys = sparse(Y);
% YYs = Ys * Ys';
% XYs = Xs * Ys';
% tmp1 = sum(sum(HH1s .* (G1s' * YYs * G1s) .* QQ1));
% tmp2 = sum(sum(HH2s .* (G2s' * YYs * G2s) .* QQ2));
% as = tmp1 + tmp2;
% 
% tmp1 = sum(sum(HH1s .* (G1s' * XYs * G1s) .* QQ1));
% tmp2 = sum(sum(HH2s .* (G2s' * XYs * G2s) .* QQ2));
% bs = 2 * (tmp1 + tmp2);
% t2 = toc(ti);
% 
% fprintf('ratio: %f\n', t1 / t2);
% equal('a', a, full(as));
% equal('b', b, full(bs));