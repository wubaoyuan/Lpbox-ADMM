%--------------------------------------------------------------------------
function params_grid_demo_clustering(task_id)
%--------------------------------------------------------------------------

base_path = '/data1/wub/lpbox_admm/';
chdir(base_path)
addpath(genpath(pwd));

% algorithm parameters
gamma_factor_list = [0.95];
rho_initial_list = [1, 1e1];
rho_change_step_list = [5, 10];
rho_increase_list = [5];
rho_upper_limit_list = [1e8];
maxIter_list = [5e3];
pNorm_list = [0.5, 1, 2, 5, 10];
data_list = [1, 2, 3, 4];

params_combination = combvec(gamma_factor_list, rho_initial_list, rho_change_step_list,...
                             rho_increase_list, rho_upper_limit_list, ...
                             maxIter_list, pNorm_list, data_list);

params_vec = params_combination(:, task_id);
gamma_factor = params_vec(1);
rho_0 = params_vec(2);
rho_change_step = params_vec(3);
mu_rho = params_vec(4);
rho_upper_limit = params_vec(5);
max_iters_ADMM = params_vec(6);
pNorm = params_vec(7);
data_index = params_vec(8);

% Lp-box ADMM parameters
params = struct('opt',2,'is_verbose', false, 'std_threshold',1e-6, ...
    'gamma_val',1.0,'gamma_factor', gamma_factor, 'initial_rho', rho_0, ...
    'x0',[],'learning_fact',1+mu_rho/100, 'rho_upper_limit', rho_upper_limit, ...
    'history_size',5, 'rel_tol',1e-6,'stop_threshold',1e-6, ...
    'max_iters',max_iters_ADMM, 'projection_lp', pNorm, 'y2_update_maxIter', inf);

if data_index == 1
   dataset_name = 'iris';
elseif data_index == 2
   dataset_name = 'wine';
elseif data_index == 3
   dataset_name = 'glass';
elseif data_index == 4
   dataset_name = 'letter20K';
end

if strcmp(dataset_name, 'letter20K')
   dataset_matrix = load('letter_only_data_20000_16.txt');
   label_gt = load('letter_label_numerical_20000_1.txt');
   params.max_iters = 3e2; 
elseif strcmp(dataset_name, 'glass') || strcmp(dataset_name, 'iris')
   dataset_matrix = load([dataset_name, '.txt']);
   label_gt = dataset_matrix(:,end);
   dataset_matrix(:,1) = []; dataset_matrix(:,end) = []; % 1st col is sample index, last col is label vector
else
   dataset_matrix = load([dataset_name, '.txt']);
   label_gt = dataset_matrix(:,1);
   dataset_matrix(:,1) = [];
end
num_cluster = length(unique(label_gt));
[num_sample, num_dimension] = size(dataset_matrix);

%% normalize each feature to [-1,1]
for i=1:num_dimension
    ma=max(dataset_matrix(:,i));
    mi=min(dataset_matrix(:,i));
    range=ma-mi;
    dataset_matrix(:,i)=2.*((dataset_matrix(:,i)-mi)./range-0.5);
end
clear range ma mi i

%% compute distance_matrix_square, based on kd-tree
% run('vl_setup')   % run this line when you firstly run this file
% matlabpool('open',16);
if strcmp(dataset_name, 'letter20K')
   num_neighbor_size = 20;   batch = 100;
else
   num_neighbor_size = 5;   batch = 1;
end
dis_matrix_kdtree_compute % return distance_matrix_square

W_matrix = log2(distance_matrix_square);
W_matrix(W_matrix==-inf) = 0;

%% - log2(distance_matrix_square), similarity matrix, degree, Laplacian matrix
M_matrix = - W_matrix;
D_matrix = diag(sum(M_matrix));
L_matrix = D_matrix - M_matrix;

%% [1]. Our ADMM Method
% min_Y  tr(Y L Y')  -  sum(diag(D_matrix)),
% s.t.      Y \in {0,1}^{K x N},
%            Y * 1_N = N/K * 1_K,  1_K^\top * Y = 1_N^\top

% L is L_matrix, K = num_cluster, N = num_sample

max_iters = 2;
record_list = zeros(max_iters, 4);
for iter = 1:max_iters
    fprintf('+ iteration [%d] ... \n',iter);
    t=tic;

    %% initialization
    if 0
	% generate random initialization
	[x0] = label_vec2binary(generate_IT_initialization(size(L_matrix,1),num_cluster),num_cluster);

    else
	% do kmeans
	%chdir('C:\Program Files\MATLAB\R2015b\toolbox\stats\stats')
	km_label_vec = kmeans(dataset_matrix,num_cluster,'start','cluster');  %run the built-in function kmeans.m
	[x0] = label_vec2binary(km_label_vec,num_cluster);
	obj_km = trace(x0'*(L_matrix)*x0)-sum(diag(D_matrix));
    end
    obj_km_list(iter) = obj_km; 
    
    params.x0 = x0;
    params.num_classes = num_cluster;
    %% solve using ADMM
    chdir('/data1/wub/lpbox_admm/functions_open/optimization')
    pwd
    [x_sol,obj_list,constraint_violation,z1,z2,time_elapsed] = ...
	ADMM_IT_clustering(L_matrix, params);
    % generate a label vector out of the binary matrix
    [~,ADMM_IT_labels] = max(x_sol,[],2);
    binary_ADMM_sol = label_vec2binary(ADMM_IT_labels,num_cluster); % the binary cluster label matrix

    %% evaluation by RAND score, compute the objective value
    [AR,RI,MI,HI] = RandIndex(label_gt,ADMM_IT_labels);
    obj_continuous = trace(x_sol'*L_matrix*x_sol)-sum(diag(D_matrix));
    obj_ADMM = trace(binary_ADMM_sol'*(L_matrix)*binary_ADMM_sol)-sum(diag(D_matrix));

    record_list(iter,:) = [RI obj_continuous obj_ADMM time_elapsed];
end

if max_iters > 1
record_list_mean = mean(record_list, 1);
record_list_std  = std(record_list, 1);
obj_km_list_mean = mean(obj_km_list, 1);
obj_km_list_std  = std(obj_km_list, 1);
else
record_list_mean = record_list;
record_list_std  = 0 .* record_list;
obj_km_list_mean = obj_km_list;
obj_km_list_std  = 0 .* obj_km_list;
end

%% save the results to file
resultStruct.record_list = record_list_mean;
resultStruct.record_list_std = record_list_std;
resultStruct.obj_km_list = obj_km_list_mean;
resultStruct.obj_km_list_std = obj_km_list_std;
resultStruct.params = params;
resultStruct.Lp = pNorm;

savePath = [base_path, '/demo/clustering/results/', dataset_name, '/'];
if ~isdir(savePath), mkdir(savePath); end
saveName = [ 'taskID_', num2str(task_id), '_', dataset_name, '_ADMM_Lp_', num2str(pNorm), '_', generate_clock_str(), '.mat' ];
save([savePath, saveName], 'resultStruct')

