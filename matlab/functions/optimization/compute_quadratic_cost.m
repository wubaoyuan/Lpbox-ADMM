%% computes the QP cost
% f(x) = x'*A*x+b'*x;

function [c] = compute_quadratic_cost(x,A,b)
c = x'*A*x+b'*x;
return;