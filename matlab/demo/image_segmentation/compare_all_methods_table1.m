% This demo implements all compared methods in Table 1 of the main script. 
% compares the different BQP methods on randomly generated graphs with random weight matrices and unary costs
% make sure that 'gco-v3.0' is compiled correctly on your machine

% Written by Baoyuan Wu and Bernard Ghanem, 2018/06/08


close all; clear all;
base_path = '/data1/wub/lpbox_admm/'; 
chdir(base_path)
addpath(genpath(pwd));

%% random graphs with random W and unary costs
node_list = 5e3; %[1e6, 5e6]; %[1e3, 5e3, 1e4, 5e4, 1e5, 5e5];  
num_trials = 5;
num_methods = 5;

% load image
im_fp = fullfile('images','cameraman.png');
I_original = im2double(imread(im_fp));
if numel(size(I_original)) == 3, I_original = rgb2gray(I_original); end;

% problem parameters
lambda = 3;

% GSA parameters
gamma = 0.5; mu = 0.1;
lg = 1.5; lm = 0.5;

% Lp-box ADMM parameters
rho_0 = 5; mu_rho = 3;   pNorm = 2; 
params = struct('opt',2,'is_verbose',false,...
    'std_threshold',1e-6,'gamma_val',1.0,'gamma_factor', 0.99,...
    'initial_rho',rho_0,'x0',[],'learning_fact',1+mu_rho/100, 'rho_upper_limit', 1000, 'history_size',5,...
    'rel_tol',1e-5,'stop_threshold',1e-3,'max_iters',1e4, 'projection_lp', pNorm);


all_results = zeros(numel(node_list),num_methods,2); % mean and std
for i = 1:length(node_list)
    n = node_list(i); 
    energy_list = zeros(num_trials,num_methods);

    % resize image + compute unary term + compute binary term
    I = imresize(I_original,sqrt(n/numel(I_original)));
    
    uParams.sig = 0.1;
    uParams.mu_b = 0.6; 
    uParams.mu_f1 = 0.2; 
    uParams.mu_f2 = uParams.mu_f1;    

    [unaryCosts, U, initLabels] = getUnaryCost(I(:), uParams); % N by 2
    unaryCosts = round(unaryCosts);
    [W,time_elapsed] = generate_binary_cost(I);
    binary_cost_mat = round(lambda*W);
    [A,b,c] = convert2Abc(unaryCosts,binary_cost_mat);

    for trial = 1:num_trials
        % setup problem
        prob_size = numel(I);
        x0 = double(rand(prob_size,1)>=0.5);
        
        % [1]. solve using GCO
        fprintf('\n\n+ processing n=%d nodes ... with GCO \n\n',prob_size);
        tic; 
        [GCO_sol,E_GCO] = compute_GCO_sol(unaryCosts,binary_cost_mat); 
        GCO_time_elapsed = toc;
        [E_GCO] = double(compute_GCO_energy(unaryCosts,binary_cost_mat,GCO_sol))
       
        % [2]. solve using penalty (GSA), the method is very slow
        fprintf('\n\n+ processing n=%d nodes ... with GSA\n\n',prob_size);
        tic; 
        [GSA_sol,fval] = GSA_BQP_unconstrained(A/2,b,gamma,lg,mu,lm,x0);
        GSA_time = toc;
        [E_GSA] = double(compute_GCO_energy(unaryCosts,binary_cost_mat,1+(GSA_sol>=0.5)));

        % [3]. solve using LP relaxation (binary constraints -> box)
        fprintf('\n\n+ processing n=%d nodes ... with LP relaxation\n\n',prob_size);
        options = optimoptions('fmincon','GradObj','on','Display','final');
        [LP_relax_sol,fval] = quadprog(A,b,[],[],[],[],zeros(prob_size,1),ones(prob_size,1),x0,options);        
        [E_LP_relax] = double(compute_GCO_energy(unaryCosts,binary_cost_mat,1+(LP_relax_sol>=0.5)));

        % [4]. solve using spectral relaxation, see citation [27] in the main manuscript
        ones_n = ones(numel(b), 1);
        A_trans_1 = A' * ones_n;
        L = [A, A_trans_1 + b; (A_trans_1 + b)', 0]./4;
        tic
        [eig_vec, eig_val] = eigs(L, 1,'sa');
        toc
        x_sol_sp = eig_vec(:, end);

        x_sol_sp(end) = [];
        x_sol_sp = (x_sol_sp + 1) ./ 2;
        E_SP = double(compute_GCO_energy(unaryCosts,binary_cost_mat,1+(x_sol_sp>=0.5)));

        % [5]. solve using lpbox ADMM
        fprintf('\n\n+ processing n=%d nodes ... with ADMM\n\n',prob_size);
        ADMM_params = params;
        ADMM_params.x0 = x0;
        [ADMM_sol,ADMM_label_vec,ADMM_obj_list,~,~,~,ADMM_time_elapsed] = ...
            ADMM_bqp_unconstrained(A/2,b,ADMM_params);  
        [E_ADMM] = double(compute_GCO_energy(unaryCosts,binary_cost_mat,1+(ADMM_sol>=0.5)))
        
        % record results
        energy_list(trial,:) = [E_GCO E_GSA E_LP_relax E_SP E_ADMM];
    end
    
    % report the mean and std energy
    if num_trials>1
        all_results(i,:,1) = mean(energy_list,1);
        all_results(i,:,2) = std(energy_list,1);
    else
        disp('got here');
        all_results(i,:,1) = energy_list;
        all_results(i,:,2) = zeros(1,num_methods);
    end
    
end

% save results
method_names = {'GCO' 'GSA' 'LP-relax' 'SP' 'ADMM'};
params.x0 = [];
if ~isdir(fullfile('./demo/image_segmentation/comparisons'))
    mkdir(fullfile('./demo/image_segmentation/comparisons'));
end 
save(fullfile('./demo/image_segmentation/comparisons',sprintf('all_results_%s.mat',generate_clock_str())),...
   'all_results','method_names','node_list','lambda','gamma','lg','mu','lm','params');


