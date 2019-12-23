%% solves a QP with box constraints using PGD
% min_x x'*A*x+b'*x   subject to: 0<=x<=1
%
% [USELESS]: Here, fA is a function handle to how A*x is computed, e.g. fA = @multiply_A

function [x,obj_list,iter,ng,alpha] = solve_QP_box_PGD(A,b,x0)
n = numel(b);
if ~exist('x0','var') || isempty(x0)
    x0 = rand(n,1);
end

std_threshold = 1e-6;
max_iters = 1e4;
iter = 1;
x = x0;
std_obj = Inf; 
obj_list = []; 
history_size = 3;
alpha = 10;

while iter<=max_iters && std_obj>=std_threshold
    % compute gradient
    grad_k = 2*(A*x)+b;
    %grad_k = 2*feval(fA,x,A)+b;
    
    % compute line search step size        
    ng = norm(grad_k);
%     alpha1 = 0.5*ng*ng/(grad_k'*A*grad_k);
    alpha = alpha*0.99;
    
%     if compute_cost(x-alpha1*grad_k,A,b)<=compute_cost(x-alpha2*grad_k,A,b)
%         alpha = alpha1;
%     else
%         alpha = alpha2;
%     end        
    
    % update x by taking a step and projecting onto the box
    x = project_box(x-alpha*grad_k);
    
    % update
    iter = iter+1;
    obj_list = [obj_list; compute_cost(x,A,b)];    
    [std_obj] = compute_std_obj(obj_list,history_size);
end
return;


%% cost function
function [c] = compute_cost(x,A,b)
c = x'*A*x+b'*x;
return;

%% computes the std of the objective history
function [std_obj] = compute_std_obj(obj_list,history_size)
std_obj = Inf;
if numel(obj_list)>=history_size
    std_obj = std(obj_list(end-history_size+1:end));
    
    % normalize 
    std_obj = std_obj/abs(obj_list(end));
end
return;

