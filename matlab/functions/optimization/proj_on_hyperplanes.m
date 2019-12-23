%% solve the problem of projecting a vector on a set of hyperplanes
% min_{x} \|x-x0\|_2^2   subject to   Ax=b

function [x,y] = proj_on_hyperplanes(A,b,x0,pcg_tol,pcg_maxit,pcg_x0)
% [1]. solve for the optimal dual variable
[y,cflag] = pcg(A*A',2*(A*x0-b),pcg_tol,pcg_maxit,[],[],pcg_x0); % could be made faster by generating AA' using symmetry

% [2]. solve for the optimal primal variable
x = x0-0.5*A'*y;
return