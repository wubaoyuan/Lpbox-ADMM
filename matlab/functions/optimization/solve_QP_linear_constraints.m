%% solve convex QP under linear constraints
% min_x x'*Ax+b'*x subject to Cx=d

function [x,y] = solve_QP_linear_constraints(A,b,C,d,x0)
n = numel(b);
neq = numel(d);

[z] = linsolve([2*A C'; C zeros(neq,neq)],[-b;d],struct('SYM',true));
%[z,cg_flag] = pcg([2*A C'; C sparse(neq,neq)],[-b;d],1e-3,1e3,[],[],[x0;zeros(neq,1)]);
%cg_flag
x = z(1:n);
y = z(n+1:end);

%norm([2*A C'; C zeros(neq,neq)]*z-[-b;d])
eig([2*A C'; C zeros(neq,neq)])
return;