%% min_x x^TAx+b^Tx such that x is {0,1}^n; Cx=d


function [best_sol,x_sol,obj_list,constraint_violation,y1,y2,time_elapsed] = ADMM_bqp_linear_eq(A,b,C,d,all_params)
%% default params
initial_params = struct('opt',1,'is_verbose',true,'stop_threshold',1e-4,'std_threshold',1e-6,'gamma_val',1.6,'gamma_factor',0.95,...
    'max_iters',1e3,'initial_rho',25,'history_size',3,'learning_fact',1+1/100,'x0',[],'pcg_tol',1e-4,'pcg_maxiters',1e3,'rel_tol',5*1e-5, 'projection_lp', 2);
if ~exist('all_params','var')
    all_params = [];
end

param_names = fieldnames(initial_params);
for p = 1:numel(param_names)    
    if ~isfield(all_params,param_names(p))
        all_params.(param_names{p}) = initial_params.(param_names{p});
    end
end


% set constants
n = numel(b);
is_verbose = all_params.is_verbose;
stop_threshold = all_params.stop_threshold;
std_threshold = all_params.std_threshold;
max_iters = all_params.max_iters;
gamma_val = all_params.gamma_val;
initial_rho = all_params.initial_rho;
learning_fact = all_params.learning_fact;  % rho increment every few iterations.
history_size = all_params.history_size; % how far back in time to see whether the cost has stabilized


%% initialization
x_sol = double(all_params.x0);
y1 = x_sol;     z1 = zeros(size(y1));
y2 = x_sol;     z2 = zeros(size(y2));
                z3 = zeros(size(d));

constraint_violation = []; % differences between x_sol, y1, and y2

% set the rho values based on Allen's paper
rho1 = initial_rho;
rho2 = rho1; 
rho3 = rho1; 

iter = 1;
obj_list = [];
[std_obj] = compute_std_obj(obj_list,history_size);
Csq = C'*C; % might not be efficient to pre-multiply depending on the size of the problem


% initiate the binary solution
prev_idx = x_sol>=0.5;
succ_change = 1;
best_sol = double(prev_idx);
best_bin_obj = compute_cost(best_sol,A,b);


%% ADMM for LBP
if is_verbose
    h = figure;
end

time_elapsed = 0;
constraint_stopping = true;

while constraint_stopping && std_obj>=std_threshold && iter<=max_iters  %&& succ_change>=all_params.rel_tol
    if is_verbose, fprintf('******* start of [iter #%d]: ****************\n',iter); end;
   
    % update y1: project onto box
    if is_verbose, fprintf('+ update y1 ... '); end;
    tic; y1 = project_box(x_sol+z1/rho1); t = toc;
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end;
    time_elapsed = time_elapsed+t;

    % update y2: project onto lp sphere
    if is_verbose, fprintf('+ update y2 ... '); end;
    tic; [y2] = project_shifted_Lp_ball(x_sol+z2/rho2,0.5*ones(n,1), all_params.projection_lp); t = toc;
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end;
    time_elapsed = time_elapsed+t;

 
    % update x: this is an exact solution to the subproblem
    % + solve a PD linear system, using pre-conditioned conjugate gradient algorithm  
    if is_verbose, fprintf('+ update x ... '); end;
    tic; 
    [x_sol,cg_flag] = pcg(rho3*Csq+sparse(1:n,1:n,rho1+rho2,n,n)+2*A,...
        -(b+z1+z2+C'*z3)+rho1*y1+rho2*y2+rho3*C'*d,...
        all_params.pcg_tol,all_params.pcg_maxiters,[],[],double(y1));
    t = toc;
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end;
    time_elapsed = time_elapsed+t;
   
 
    % update z1 and z2 and z3
    if is_verbose, fprintf('+ update z1 and z2 and z3 ... '); end;
    tic; 
    z1 = z1+gamma_val*rho1*(x_sol-y1);
    z2 = z2+gamma_val*rho2*(x_sol-y2); 
    z3 = z3+gamma_val*rho3*(C*x_sol-d);
    t = toc;
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end;
    time_elapsed = time_elapsed+t;
   
    
    % evaluate this iteration
    temp_res = [norm(x_sol-y1) norm(x_sol-y2)]/n;
    constraint_stopping = (max(temp_res)>=stop_threshold);
    constraint_violation = [constraint_violation; temp_res];     

    % try a different stopping criterion
    cur_idx = x_sol>=0.5;
    if iter> 10; succ_change = sum(cur_idx~=prev_idx)/n; end;
    prev_idx = cur_idx;
    cur_obj = compute_cost(prev_idx,A,b);
    obj_list = [obj_list; compute_cost(x_sol,A,b)];
    [std_obj] = compute_std_obj(obj_list,history_size);

    
    if is_verbose, fprintf('******* end of [iter #%d]: std = [%3.3f]; time = [%3.3f] ****************\n',iter,log10(std_obj),time_elapsed); end;
    iter = iter + 1;    
    
    % increase rhos and update gamma is needed
    if mod(iter,5)==0
        rho1 = learning_fact*rho1;
        rho2 = learning_fact*rho2;
        rho3 = learning_fact*rho3;
        gamma_val = max(gamma_val*all_params.gamma_factor,1);        
    end
    
    % maintain best binary solution so far; in case the cost function oscillates
    if best_bin_obj > cur_obj && max(temp_res)>=1e-3
        best_bin_obj = cur_obj;
        best_sol = x_sol;
    end
   
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
    std_obj = std_obj/(max(abs(obj_list(end)),eps));
end
return;
