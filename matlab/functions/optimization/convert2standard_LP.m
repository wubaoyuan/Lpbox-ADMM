%% converts an ordinary LP formulated as: 
% min_x a^Tx subject to Bx=c; Dx<=e
%  
% to standard form:
% min_z f^Tz subject to Gz=h; z>=0
%
% This conversion is possible by adding auxiliary variables:
% x = xp-xm; s=e-Dx; z = [xp; xm; s]

function [f,G,h] = convert2standard_LP(a,B,c,D,e)
% covnert the variables by adding axiliary variables
neq = size(B,1);
nineq = numel(e);

f = [a; -a];
G = [B -B zeros(neq,nineq); D -D eye(nineq)];
h = [c;e];
return