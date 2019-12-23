%% min_x (x^T A x + b^T x) such that x is {0,1}^n
% ADMM update steps with x0 being feasible and binary
% try to use sparse matrix for A

function [best_sol,x_sol,obj_list,constraint_violation,y1,y2,time_elapsed] = ADMM_bqp_unconstrained(A,b,all_params)
%% default params
initial_params = struct('opt',2,'is_verbose',true,'stop_threshold',1e-4,'std_threshold',1e-5,'gamma_val',1.6,'gamma_factor',0.95,...
    'max_iters',1e3,'initial_rho',25,'history_size',3,'learning_fact',1+1/100,'x0',[],'pcg_tol',1e-3,'pcg_maxiters',1e3,...
    'rel_tol',5*1e-5,'imsize',[],'save_dir',[], 'projection_lp', 2, 'z2_update_maxIter', inf, 'rho_change_step', 5);

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
opt = all_params.opt;
is_verbose = all_params.is_verbose;
stop_threshold = all_params.stop_threshold;
std_threshold = all_params.std_threshold;
max_iters = all_params.max_iters;
initial_rho = all_params.initial_rho;
rho_change_step = all_params.rho_change_step;
gamma_val = all_params.gamma_val;
learning_fact = all_params.learning_fact;  % rho increment every few iterations. 
history_size = all_params.history_size; % how far back in time to see whether the cost has stabilized

%% initialization
x_sol = single(all_params.x0);
y1 = x_sol; z1 = zeros(size(y1),'single');
y2 = x_sol; z2 = zeros(size(y2),'single');
constraint_violation = []; % differences between x_sol, y1, and y2
rho1 = initial_rho; 
rho2 = rho1;
obj_list = [];
[std_obj] = compute_std_obj(obj_list,history_size);

% invert the matrix assuming constant rho
if opt == 1
    if is_verbose, fprintf('+ computing inverse of the core matrix ... '); end
    tic; invAext = 0.5*inv(A+rho*eye(n)); t = toc;
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end
end

% initiate the binary solution
prev_idx = x_sol>=0.5;
succ_change = 1;
best_sol = double(prev_idx);
best_bin_obj = compute_cost(best_sol,A,b);


%% ADMM 
if ~exist(all_params.save_dir,'dir') && ~isempty(all_params.imsize)
    if ~exist(fullfile(all_params.save_dir,'x'),'dir')
        mkdir(all_params.save_dir,'x');
    end
    if ~exist(fullfile(all_params.save_dir,'y1'),'dir')
        mkdir(all_params.save_dir,'y1');
    end
    if ~exist(fullfile(all_params.save_dir,'y2'),'dir')
        mkdir(all_params.save_dir,'y2');
    end
    
    imwrite(reshape(x_sol,all_params.imsize),...
        fullfile(all_params.save_dir,'x',sprintf('%03d_x_sol.png',0)));
    imwrite(reshape(y1,all_params.imsize),...
        fullfile(all_params.save_dir,'y1',sprintf('%03d_y1.png',0)));
    imwrite(reshape(y2,all_params.imsize),...
        fullfile(all_params.save_dir,'y2',sprintf('%03d_y2.png',0)));
end

if is_verbose 
    h = figure;
end

time_elapsed = 0;
constraint_stopping = true;
iter = 1; 
while constraint_stopping && std_obj>=std_threshold && iter<=max_iters %&& succ_change>=all_params.rel_tol
    if is_verbose, fprintf('******* start of [iter #%d]: ****************\n',iter); end;

    % update y1
    if is_verbose, fprintf('+ update y1 ... '); end;
    tic; y1 = project_box(x_sol+z1/rho1); t = toc;
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end;
    time_elapsed = time_elapsed+t;

    % update y2
    if is_verbose, fprintf('+ update y2 ... '); end;
    tic; [y2] = project_shifted_Lp_ball(x_sol+z2/rho2,0.5*ones(n,1), all_params.projection_lp); t = toc;
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end;
    time_elapsed = time_elapsed+t;  
 
    % update x: this is an exact solution to the subproblem
    if is_verbose, fprintf('+ update x ... '); end;
    tic; rhs = rho1*y1+rho2*y2-(b+z1+z2);

    switch opt
        case 1 % using matrix inversion
            x_sol = invAext*rhs; 
            
        case 2 % using PCG, DEFAULT
            [x_sol,cg_flag] = pcg((2*A+sparse(1:n,1:n,rho1+rho2,n,n)),double(rhs),all_params.pcg_tol,all_params.pcg_maxiters,[],[],double(y1));
    end
    t = toc;
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end;
    time_elapsed = time_elapsed+t;
    
    % update z1 and z2
    if is_verbose, fprintf('+ update z1 and z2 ... '); end;
    tic;     
    z1 = z1+gamma_val*rho1*(x_sol-y1); 
    z2 = z2+gamma_val*rho2*(x_sol-y2); 
    t = toc; 
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end;
    time_elapsed = time_elapsed+t;
   
    
    % evaluate this iteration
    temp_res = [norm(x_sol-y1) norm(x_sol-y2)]/max(norm(x_sol),eps);
    constraint_stopping = max(temp_res)>=stop_threshold;
    constraint_violation = [constraint_violation; temp_res];
    obj_list = [obj_list; compute_cost(x_sol,A,b)];
    [std_obj] = compute_std_obj(obj_list,history_size);
    
    % try a different stopping criterion
    cur_idx = x_sol>=0.5;
    if iter> 10; succ_change = sum(cur_idx~=prev_idx)/numel(cur_idx); end;
    prev_idx = cur_idx;
    cur_obj = compute_cost(prev_idx,A,b);
    
    % maintain best binary solution so far; in case the cost function oscillates
    if best_bin_obj >= cur_obj
        best_bin_obj = cur_obj;
        best_sol = x_sol;
    end
    
    if is_verbose, fprintf('******* end of [iter #%d]: std = [%3.2f]; time = [%3.1f]; relc = [%3.1f]; sc=[%3.1f] ****************\n', ...
       iter,log10(std_obj),time_elapsed,log10(succ_change),log10(max(temp_res))); end;    
    
    % increase rho1 and rho2
    if mod(iter, rho_change_step)==0
        rho1 = learning_fact*rho1;
        rho2 = learning_fact*rho2;
        gamma_val = max(gamma_val*all_params.gamma_factor,1);
    end  
    
    % save images if needed
    if exist(all_params.save_dir,'dir') && ~isempty(all_params.imsize)        
        imwrite(reshape(x_sol,all_params.imsize),...
            fullfile(all_params.save_dir,'x',sprintf('%03d_x_sol.png',iter)));
        imwrite(reshape(y1,all_params.imsize),...
            fullfile(all_params.save_dir,'y1',sprintf('%03d_y1.png',iter)));
        imwrite(reshape(y2,all_params.imsize),...
            fullfile(all_params.save_dir,'y2',sprintf('%03d_y2.png',iter)));
    end
    
    % increment counter
    iter = iter + 1;
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
