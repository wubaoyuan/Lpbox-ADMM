%% min_x  x'*A*x+b'*x+gamma*x'*(1-x)-mu*sum_i (ln(x_i)+ln(1-x_i))  s.t. 0<=x<=1 
%
% [1]. Walter Murray and Kien-Ming Ng, "An algorithm for nonlinear
% optimization problems with binary variables", Computer Optim Appl (2010)


function [x_sol,fval] = GSA_BQP_unconstrained(A,b,gamma,lg,mu,lm,x0)
n = numel(b);
obj_list = [];rel_list = []; do_stop = false;
options = optimoptions('fmincon','GradObj','on','Display','final',...
    'TolFun',1e-4,'TolX',1e-5);

while ~do_stop
    [x_sol,fval] = fmincon(@(x)cost(x,A,b,gamma,mu),x0,[],[],[],[],zeros(n,1),ones(n,1),[],options);
    obj_list = [obj_list; x_sol'*A*x_sol+b'*x_sol];
    rel_list = [rel_list; norm(x_sol-x0)];    
    x0 = x_sol;
    
    % change the tradeoff parameters
    gamma = lg*gamma;
    mu = lm*mu;
    
    % stopping criterion
    do_stop = gamma>=100;
end

return;

function [c,g] = cost(x,A,b,gamma,mu)
c = x'*A*x+b'*x+gamma*(sum(x)-x'*x)-mu*(sum(log(x)+log(1-x)));
g = 2*A*x+b+gamma-2*gamma*x-mu*(1./x+1./(1-x));
return;