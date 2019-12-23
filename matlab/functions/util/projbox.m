function [x] = projbox(x)
x(x<0)=0;
x(x>1)=1;