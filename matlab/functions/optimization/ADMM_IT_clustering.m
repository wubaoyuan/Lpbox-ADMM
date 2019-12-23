%% min_x \sum_{i=1}^c (x_i^T A x_i) such that 
% + x is {0,1}^{n * c}; 
% + sum_i x_i=1 (each instance is assigned to one and only one cluster)
% + sum_j x_ij=n/c (each cluster has the same number of members)
% 
% Here, 
% + c is the number of clusters; 
% + A is the Laplacian of the graph
% + x0 is an initial guess (in matrix form: nxc values)
% + ADMM update steps with x0 being feasible and binary;
% + try to use sparse matrix for A.


function [x_sol,obj_list,constraint_violation,y1,y2,time_elapsed] = ADMM_IT_clustering(A,all_params)
%% default params
startt = tic;
initial_params = struct('opt',1,'num_classes',2,'is_verbose',true,'stop_threshold',1e-3,'std_threshold',1e-5,'gamma_val',1.6,'gamma_factor',0.95,...
    'max_iters',1e3,'initial_rho',25, 'rho_upper_limit', inf, 'history_size',3,'learning_fact',1+1/100,'x0',[], ...
    'pcg_tol',1e-3,'pcg_maxiters',1e3,'rel_tol',5*1e-5,'imsize',[], 'projection_lp', 2, 'z2_update_maxIter', inf);
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
num_nodes   = size(A,1);
num_classes = all_params.num_classes;
opt = all_params.opt;
is_verbose = all_params.is_verbose;
stop_threshold = all_params.stop_threshold;
std_threshold = all_params.std_threshold;
max_iters = all_params.max_iters;
initial_rho = all_params.initial_rho;
learning_fact = all_params.learning_fact;  % rho increment every few iterations. This part needs to be replaced with a less heuristic and more systematic update operator
history_size = all_params.history_size; % how far back in time to see whether the cost has stabilized
gamma_val = all_params.gamma_val;

%% initialization
if isempty(all_params.x0)
    all_params.x0 = zeros(num_nodes,num_classes);
end

x_sol = double(all_params.x0(:));
y1 = x_sol; z1 = zeros(size(y1));
y2 = x_sol; z2 = zeros(size(y2));
            z3 = zeros(num_nodes,1);
            z4 = zeros(num_classes,1);

constraint_violation = Inf(2,1); % differences between x_sol, y1, and y2

rho1 = initial_rho; 
rho2 = rho1;
rho3 = rho1;
rho4 = rho1;

iter = 1;
obj_list = [];
obj_list(1) = compute_cost(x_sol,A,num_classes);
[std_obj] = compute_std_obj(obj_list,history_size);
time_elapsed = 0;

% current index
[~,prev_idx] = max(reshape(x_sol,[num_nodes num_classes]),[],2);
succ_change = 1;


%% ADMM loop
ideal_cluster_size = num_nodes/num_classes;
constraint_stopping = true;

if is_verbose
    h = figure;
end

while constraint_stopping && iter<=max_iters && std_obj>=std_threshold %&& succ_change>=all_params.rel_tol
    if is_verbose, fprintf('******* start of [iter #%d]: ****************\n',iter); end;
   
    %% [1]. update y1 (try over-relaxation?)
    if is_verbose, fprintf('+ update y1 ... '); end;
    tic; y1 = project_box(x_sol+z1/rho1); t = toc;
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end;
    time_elapsed = time_elapsed+t;

    %% [1]. update y2 (try over-relaxation?)
    if is_verbose, fprintf('+ update y2 ... '); end;
    tic; [y2] = project_shifted_Lp_ball(x_sol+z2/rho2,0.5*ones(num_classes*num_nodes,1), all_params.projection_lp); t = toc;
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end;
    time_elapsed = time_elapsed+t;

 
    %% [2]. update x: this is an exact solution to the subproblem (try approximating it)
    % + solve a large and sparse convex QP with equality constraints
    if is_verbose, fprintf('+ update x ... '); end;
    
    switch opt            
        case 1 % using matrix inversion
            x_sol = invAext*rhs;
    
        case 2 % separate the summation equality and the cluster uniformity constraint from the rest
            tic; 
            rhs = (rho1*y1+rho2*y2)-(z1+z2)+rho3+rho4*ideal_cluster_size;
            for p = 1:num_classes
                rhs(1+(p-1)*num_nodes:p*num_nodes) = rhs(1+(p-1)*num_nodes:p*num_nodes)-z3-z4(p);
            end
            
            [x_sol,cg_flag,~,pcg_iter] = pcg(@afun2,rhs,all_params.pcg_tol,all_params.pcg_maxiters,[],[],...
                double(y1),A,num_nodes,num_classes,rho1,rho2,rho3,rho4);            
    end
           
        
    t = toc;
    
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end;
    time_elapsed = time_elapsed+t;
    
    %% [3]. update z1 - z4 (try over-relaxation?)
    reshape_x = reshape(x_sol,[num_nodes num_classes]);
    
    if is_verbose, fprintf('+ update z1 to z4 ... '); end;
    tic; 
    z1 = z1+gamma_val*rho1*(x_sol-y1);
    if iter < all_params.z2_update_maxIter
        z2 = z2+gamma_val*rho2*(x_sol-y2); 
    else
        z2 = z2;
    end
    z3 = z3+gamma_val*rho3*(sum(reshape_x,2)-1);
    z4 = z4+gamma_val*rho4*vec(sum(reshape_x,1)-ideal_cluster_size);
    t = toc;
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end;
    time_elapsed = time_elapsed+t;
    
    
    %% evaluate this iteration
    constraint_violation = [norm(x_sol-y1) norm(x_sol-y2)]/max(norm(x_sol),eps);
    obj_list = [obj_list; compute_cost(x_sol,A,num_classes)];
    [std_obj] = compute_std_obj(obj_list,history_size);
    
    % try a different stopping criterion
    [~,cur_idx] = max(reshape(x_sol,[num_nodes num_classes]),[],2);
    succ_change = sum(cur_idx~=prev_idx)/numel(cur_idx);
    prev_idx = cur_idx;
    
    if is_verbose, fprintf('******* end of [iter #%d]: std = [%3.2f]; time = [%3.1f]; relc = [%3.1f] ****************\n',iter,log10(std_obj),time_elapsed,log10(succ_change)); end;
    iter = iter + 1;    
    
    % increase rho1 and rho2    
    if mod(iter,5)==0    %&& iter < all_params.z2_update_maxIter
        if rho1 <= all_params.rho_upper_limit
            rho1 = learning_fact*rho1;
        end
        if rho2 <= all_params.rho_upper_limit
            rho2 = learning_fact*rho2;
        end
        if rho3 <= all_params.rho_upper_limit
             if opt==3; rho3 = learning_fact*rho3; end;
        end
        
        % make gamma smaller for more stability close to the solution
        gamma_val = max(gamma_val*all_params.gamma_factor,1);
        
        % plot some results
        if is_verbose
            figure(h); 
            if ~isempty(all_params.imsize)
                imshow(reshape(cur_idx,all_params.imsize),[]);
                title(sprintf('obj=[%3.3f]',log10(obj_list(end))));
            else
                ncs = min(num_classes,6);
                for c = 1:ncs
                    subplot(1,ncs,c); plot(sort(x_sol(1+(c-1)*num_nodes:c*num_nodes)),'.'); 
                    ylim([-0.1 1.1]); xlim([1 num_nodes]); title(sprintf('[c%d]. obj=[%3.3f]',c,log10(obj_list(end))));
                end
            end
        end
    end
end

% reshape the solution as the size of x0
x_sol = reshape(x_sol,size(all_params.x0));
y1 = reshape(y1,size(all_params.x0));
y2 = reshape(y2,size(all_params.x0));
time_elapsed = toc(startt);
return;


%% function to multiply a vector with the QP matrix
function [res] = afun2(x,A,num_nodes,num_classes,rho1,rho2,rho3,rho4)
reshape_x = reshape(x,[num_nodes num_classes]);
res = bsxfun(@plus,2*A*reshape_x+(rho1+rho2)*reshape_x,rho3*sum(reshape_x,2));

% add to each column a different vector
temp_sum = sum(reshape_x,1);
for p = 1:num_classes
    res(:,p) = res(:,p) + rho4*temp_sum(p);
end
res = res(:);
return;


%% cost function (binary)
function [c] = compute_cost(x,A,num_classes)
num_nodes = size(A,1); 
c = 0;
% x = double(x>=0.5);
for p = 1:num_classes
    c = c + x(1+(p-1)*num_nodes:p*num_nodes)'*(A*x(1+(p-1)*num_nodes:p*num_nodes));
end
% c = c+x'*b;
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
