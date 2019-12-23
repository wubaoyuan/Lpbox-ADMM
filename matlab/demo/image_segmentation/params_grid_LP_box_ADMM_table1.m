%--------------------------------------------------------------------------
function params_grid_LP_box_ADMM_table1(task_id)
%--------------------------------------------------------------------------

base_path = '/data1/wub/lpbox_admm/'; 
chdir(base_path)
addpath(genpath(pwd));

chdir('./demo/image_segmentation')


% load image
im_fp = fullfile('images','cameraman.png');
I_original = im2double(imread(im_fp));
if numel(size(I_original)) == 3, I_original = rgb2gray(I_original); end;

% problem parameters
lambda = 3;

% algorithm parameters 
gamma_factor_list = [0.95, 0.99];
rho_initial_list = [1e-2, 1e-1, 1, 1e1, 5e1];
rho_change_step_list = [5, 10]; 
rho_increase_list = [1, 5, 10, 20];
rho_upper_limit_list = [1e4, 1e6]; 
maxIter_list = 10;%[1.5e4]; 
pNorm_list = 10; %[0.5, 1, 2, 5, 10]; 
node_list = 5e3; %[5e3, 1e4, 5e4, 1e5, 5e5, 1e6, 5e6];  

params_combination = combvec(gamma_factor_list, rho_initial_list, rho_change_step_list, ...
                             rho_increase_list, rho_upper_limit_list, ...
                             maxIter_list, pNorm_list, node_list); 
                                        
params_vec = params_combination(:, task_id);	  
gamma_factor = params_vec(1); 
rho_0 = params_vec(2); 
rho_change_step = params_vec(3); 
mu_rho = params_vec(4); 
rho_upper_limit = params_vec(5); 
max_iters = params_vec(6); 
pNorm = params_vec(7); 
num_nodes = params_vec(8); 
						  
% Lp-box ADMM parameters
params = struct('opt',2,'is_verbose',true, 'std_threshold',1e-6, ...
    'gamma_val',1.0,'gamma_factor', gamma_factor, 'initial_rho', rho_0, ...
    'x0',[],'learning_fact',1+mu_rho/100, 'rho_upper_limit', rho_upper_limit, ...
    'history_size',5, 'rel_tol',1e-5,'stop_threshold',1e-3, ...
    'max_iters',max_iters, 'projection_lp', pNorm, 'y2_update_maxIter', inf);

% resize image + compute unary term + compute binary term
I = imresize(I_original,sqrt(num_nodes/numel(I_original)));

uParams.sig = 0.1;
uParams.mu_b = 0.6; 
uParams.mu_f1 = 0.2; 
uParams.mu_f2 = uParams.mu_f1;    

[unaryCosts, U, initLabels] = getUnaryCost(I(:), uParams); % N by 2
unaryCosts = round(unaryCosts);
[W,time_elapsed] = generate_binary_cost(I);
binary_cost_mat = round(lambda*W);
[A,b,c] = convert2Abc(unaryCosts,binary_cost_mat);

% setup problem
prob_size = numel(I);
x0 = double(rand(prob_size,1)>=0.5);

% [4]. solve using ADMM
fprintf('\n\n+ processing n=%d nodes ... with ADMM\n\n',prob_size);
ADMM_params = params;
ADMM_params.x0 = x0;
[ADMM_sol,ADMM_label_vec,ADMM_obj_list,~,~,~,ADMM_time_elapsed] = ...
    ADMM_bqp_unconstrained(A/2,b,ADMM_params);  
[E_ADMM] = double(compute_GCO_energy(unaryCosts,binary_cost_mat,1+(ADMM_sol>=0.5)));

result_struct_admm = struct('energy', E_ADMM, ...  %'best_solution', ADMM_sol, ...
                            'final_solution', ADMM_label_vec, ...
                            'obj_curve', ADMM_obj_list, ...
                            'runtime', ADMM_time_elapsed, ...
                            'params', ADMM_params);

%--------------------------- save results                    
save_name = ['taskID_', num2str(task_id), '_result_lpbox_admm_n_', num2str(num_nodes), ...
             '_pNorm_', num2str(pNorm), '_', generate_clock_str(), '.mat'];

if ~isdir(fullfile(base_path, 'demo/image_segmentation', 'result'))
  mkdir(fullfile(base_path, 'demo/image_segmentation', 'result')); 
end
save_path = fullfile(base_path, 'demo/image_segmentation', 'result');

save(fullfile(save_path, save_name));
