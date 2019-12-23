%% min_x \sum_{i=1}^c (x_i^T A x_i) + b_i^T x_i) such that x is {0,1}^nc
% 
% Here, 
% + c is the number of classes; 
% + A is the Laplacian of the graph
% + b is the unary vector for each class concatenated 
% + x0 is an initial guess (in matrix form: nxc values)
% + ADMM update steps with x0 being feasible and binary;
% + try to use sparse matrix for A.


function [x_sol,obj_list,constraint_violation,z1,z2,time_elapsed] = ADMM_bqp_exact_multiclass(A,b,all_params)
%% default params
initial_params = struct('opt',1,'is_verbose',true,'stop_threshold',1e-3,'std_threshold',1e-5,'gamma_val',1.6,'gamma_factor',0.95,...
    'max_iters',1e3,'initial_rho',25,'history_size',3,'learning_fact',1+1/100,'x0',[],'pcg_tol',1e-3,'pcg_maxiters',1e3,...
    'rel_tol',5*1e-5,'imsize',[],'save_dir',[]);
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
[num_nodes,num_classes] = size(all_params.x0);
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
x_sol = double(all_params.x0(:));
z1 = x_sol; y1 = zeros(size(z1));
z2 = x_sol; y2 = zeros(size(z2));
constraint_violation = Inf(2,1); % differences between x_sol, z1, and z2
rho1 = initial_rho; 
rho2 = rho1;
iter = 1;
obj_list = [];
[std_obj] = compute_std_obj(obj_list,history_size);
time_elapsed = 0;
%pcg_params = struct('tol',1e-3,'maxit',1e2);

if opt == 3
    y3 = zeros(num_nodes,1);
    rho3 = rho1;
end

% current index
[~,prev_idx] = max(reshape(x_sol,[num_nodes num_classes]),[],2);
succ_change = 1;

%% ADMM loop
constraint_stopping = true;

if ~exist(all_params.save_dir,'dir') && ~isempty(all_params.imsize)
    if ~exist(fullfile(all_params.save_dir,'x'),'dir')
        mkdir(all_params.save_dir,'x');
    end
    if ~exist(fullfile(all_params.save_dir,'z1'),'dir')
        mkdir(all_params.save_dir,'z1');
    end
    if ~exist(fullfile(all_params.save_dir,'z2'),'dir')
        mkdir(all_params.save_dir,'z2');
    end
    
    % save the initial solutions
    %{
    [~,temp_idx] = max(reshape(x_sol,[num_nodes num_classes]),[],2);
    imwrite(reshape(temp_idx,all_params.imsize),gray(num_classes),...
            fullfile(all_params.save_dir,'x',sprintf('%03d_x_sol.png',0)));    
    
    [~,temp_idx] = max(reshape(z1,[num_nodes num_classes]),[],2);
    imwrite(reshape(temp_idx,all_params.imsize),gray(num_classes),...
            fullfile(all_params.save_dir,'z1',sprintf('%03d_z1.png',0)));    

    [~,temp_idx] = max(reshape(z2,[num_nodes num_classes]),[],2);
    imwrite(reshape(temp_idx,all_params.imsize),gray(num_classes),...
            fullfile(all_params.save_dir,'z2',sprintf('%03d_z2.png',0)));
    %}
end

if is_verbose
    h = figure;
end

%while constraint_stopping && iter<=max_iters && std_obj>=std_threshold
while constraint_stopping && iter<=max_iters && std_obj>=std_threshold && succ_change>=all_params.rel_tol
    if is_verbose, fprintf('******* start of [iter #%d]: ****************\n',iter); end;
    
    %% [1]. update x: this is an exact solution to the subproblem (try approximating it)
    % + solve a large and sparse convex QP with equality constraints
    if is_verbose, fprintf('+ update x ... '); end;
    
    switch opt
        case 1          
            tic; rhs = [(rho1*z1+rho2*z2)-(b+y1+y2); ones(num_nodes,1)];    
            [qp_sol,cg_flag,~,pcg_iter] = pcg(@afun,rhs,all_params.pcg_tol,all_params.pcg_maxiters,[],[],rand(num_nodes*(num_classes+1),1),A,num_nodes,num_classes,rho1,rho2);
            %[qp_sol,cg_flag] = pcg(@(x)afun(x,A),rhs,pcg_params.tol,pcg_params.maxit,[],[],[],A);
            x_sol = qp_sol(1:num_classes*num_nodes);

            cg_flag
            pcg_iter
        case 2
            tic; 
            total_rhs = (rho1*z1+rho2*z2)-(b+y1+y2);
            for p = 1:num_classes
                % solve the graph cut problem for each class
                [x_sol(1+(p-1)*num_nodes:p*num_nodes),cg_flag] = pcg((2*A+sparse(1:num_nodes,1:num_nodes,rho1+rho2,num_nodes,num_nodes)),...
                    total_rhs(1+(p-1)*num_nodes:p*num_nodes),all_params.pcg_tol,all_params.pcg_maxiters,[],[],z1(1+(p-1)*num_nodes:p*num_nodes));
                cg_flag
            end

            % make sure the solutions for each node sum to 1
            reshape_x = reshape(x_sol,[num_nodes,num_classes]);
            x_sol = bsxfun(@times,reshape_x,min(eps,1./sum(reshape_x,2)));
            x_sol = x_sol(:);   
            
        case 3 % separate the summation equality from the rest
            tic; 
            rhs = (rho1*z1+rho2*z2)-(b+y1+y2)+rho3;
            for p = 1:num_classes
                rhs(1+(p-1)*num_nodes:p*num_nodes) = rhs(1+(p-1)*num_nodes:p*num_nodes)-y3;
            end
            
            [x_sol,cg_flag,~,pcg_iter] = pcg(@afun2,rhs,all_params.pcg_tol,all_params.pcg_maxiters,[],[],...
                double(z1),A,num_nodes,num_classes,rho1,rho2,rho3);  % initiliaze with z1 or z2 instead of x_sol        
    end
           
        
    t = toc;
    
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end;
    time_elapsed = time_elapsed+t;
    
    %% [2]. update z1 (try over-relaxation?)
    if is_verbose, fprintf('+ update z1 ... '); end;
    tic; z1 = project_box(x_sol+y1/rho1); t = toc;
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end;
    time_elapsed = time_elapsed+t;
    
    %% [2]. update z2 (try over-relaxation?)
    if is_verbose, fprintf('+ update z2 ... '); end;
    tic; [z2] = project_shifted_L2_ball(x_sol+y2/rho2,0.5*ones(num_classes*num_nodes,1)); t = toc;
    
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end;
    time_elapsed = time_elapsed+t;
    
    %% [3]. update y1 and y2 (try over-relaxation?)
    if is_verbose, fprintf('+ update y1 and y2 ... '); end;
    tic; y1 = y1+gamma_val*rho1*(x_sol-z1);
    y2 = y2+gamma_val*rho2*(x_sol-z2); t = toc;
    if is_verbose, fprintf('finished in [%3.3f] seconds \n',t); end;
    time_elapsed = time_elapsed+t;
    
    if opt==3
        % update y3
        y3 = y3+gamma_val*rho3*(sum(reshape(x_sol,[num_nodes num_classes]),2)-1);
    end
    
    %% evaluate this iteration
    constraint_violation = [norm(x_sol-z1) norm(x_sol-z2)]/max(norm(x_sol),eps);
    constraint_stopping = max(constraint_violation)>=stop_threshold;
    obj_list = [obj_list; compute_cost(x_sol,A,b)];
    [std_obj] = compute_std_obj(obj_list,history_size);
    
    % try a different stopping criterion
    [~,cur_idx] = max(reshape(x_sol,[num_nodes num_classes]),[],2);
    succ_change = sum(cur_idx~=prev_idx)/numel(cur_idx);
    prev_idx = cur_idx;
    
    if is_verbose, fprintf('******* end of [iter #%d]: std = [%3.2f]; time = [%3.1f]; relc = [%3.1f] ****************\n',iter,log10(std_obj),time_elapsed,log10(succ_change)); end;    
    
    % increase rho1 and rho2    
    if mod(iter,5)==0
        rho1 = learning_fact*rho1;
        rho2 = learning_fact*rho2;
        if opt==3; rho3 = learning_fact*rho3; end;
        
        % make gamma smaller for more stability close to the solution
        gamma_val = max(gamma_val*all_params.gamma_factor,1);
        
        % plot some results
        if is_verbose
            figure(h); 
            if ~isempty(all_params.imsize)
                imshow(reshape(cur_idx,all_params.imsize),[]);
                title(sprintf('obj=[%3.3f]',log10(obj_list(end))));
            else
                for c = 1:num_classes
                    subplot(1,num_classes,c); plot(sort(x_sol(1+(c-1)*num_nodes:c*num_nodes)),'.'); 
                    ylim([-0.1 1.1]); xlim([1 num_nodes]); title(sprintf('obj=[%3.3f]',log10(obj_list(end))));
                end
            end
        end
    end
    
    % save images if needed
    if exist(all_params.save_dir,'dir') && ~isempty(all_params.imsize)  
        [~,temp_idx] = max(reshape(x_sol,[num_nodes num_classes]),[],2);
        imwrite(reshape(temp_idx,all_params.imsize),gray(num_classes),...
                fullfile(all_params.save_dir,'x',sprintf('%03d_x_sol.png',iter)));    

        [~,temp_idx] = max(reshape(z1,[num_nodes num_classes]),[],2);
        imwrite(reshape(temp_idx,all_params.imsize),gray(num_classes),...
                fullfile(all_params.save_dir,'z1',sprintf('%03d_z1.png',iter)));    

        [~,temp_idx] = max(reshape(z2,[num_nodes num_classes]),[],2);
        imwrite(reshape(temp_idx,all_params.imsize),gray(num_classes),...
                fullfile(all_params.save_dir,'z2',sprintf('%03d_z2.png',iter)));
    end  
    
    % increment counter
    iter = iter + 1;
end

% reshape the solution as the size of x0
x_sol = reshape(x_sol,size(all_params.x0));
return;


%% function to multiply a vector with the QP matrix
function [res] = afun(x,A,num_nodes,num_classes,rho1,rho2)
reshape_x = reshape(x(1:num_nodes*num_classes),[num_nodes num_classes]);
res = bsxfun(@plus,2*A*reshape_x+(rho1+rho2)*reshape_x,x(num_nodes*num_classes+1:end));
res = [res(:); sum(reshape_x,2)];
return;

function [res] = afun2(x,A,num_nodes,num_classes,rho1,rho2,rho3)
reshape_x = reshape(x,[num_nodes num_classes]);
res = bsxfun(@plus,2*A*reshape_x+(rho1+rho2)*reshape_x,rho3*sum(reshape_x,2));
res = res(:);
return;


% function [res] = afun2(x,A)
% num_nodes = size(A,1);
% num_classes = numel(x)/num_nodes;
% reshape_x = reshape(x(1:num_nodes*num_classes),[num_nodes num_classes]);
% res = 2*A*reshape_x;
% res = res(:);
% return;


%% cost function (binary)
function [c] = compute_cost(x,A,b)
num_nodes = size(A,1); num_classes = numel(b)/num_nodes;
bterm = 0;
% x = double(x>=0.5);
for p = 1:num_classes
    bterm = bterm + x(1+(p-1)*num_nodes:p*num_nodes)'*(A*x(1+(p-1)*num_nodes:p*num_nodes));
end
c = b'*x+bterm;
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