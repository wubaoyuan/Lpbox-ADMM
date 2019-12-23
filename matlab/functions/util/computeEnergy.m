function [fobj] = computeEnergy(x,A,b,c)
x=x(:);
x =  round(x);
fobj = 0.5*x'*A*x+b'*x+c;

    