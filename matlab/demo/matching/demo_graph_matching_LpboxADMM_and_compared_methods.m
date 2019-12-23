%{
     A demo comparison of different graph matching methods on the on CMU House dataset. It will reproduce the results reported in Figure 6 in our pami manuscript.
     The original file is written by Feng Zhou (zhfe99@gmail.com). We modify it by adding our Lp-box ADMM algorithm, to compare with different matching methods.
     -- Baoyuan Wu, Bernard Ghanem,  2018/06/08
%}

clear; clear all; 
chdir('/data1/wub/lpbox_admm/')
addpath(genpath(pwd))

% clear variables;
prSet(1);
tag = 'house';
num_frames = 111; % there are 111 frames in the house data (video)

% src parameter
%% randomly remove 'num_node_remove' nodes ( [30 - num_node_remove] nodes left ) in both images, then not every one node has a matching node
num_node_remove = 5; 
nIn = [30 30] - num_node_remove; 

frame_gap_list = [1, 10:10:90];
for frame_gap_id = 1: length(frame_gap_list)
    frame_gap = frame_gap_list(frame_gap_id);
    num_pairs = num_frames-frame_gap; 

    %% define the results matrix of all methods
    % one row is one pair, 1st column records the accuracy, 2nd column saves objective value, 3rd column saves the running time
    results_matrix_PM = zeros(num_pairs, 3); 
    results_matrix_SM = zeros(num_pairs, 3);
    results_matrix_IPFP_U = zeros(num_pairs, 3);
    results_matrix_IPFP_S = zeros(num_pairs, 3);
    results_matrix_RRWM = zeros(num_pairs, 3);
    results_matrix_FGM_U = zeros(num_pairs, 3);
    results_matrix_FGM_D = zeros(num_pairs, 3);
    results_matrix_QP = zeros(num_pairs, 3);
    
    pNormList = [0.5, 1, 2, 5, 10]; 
    numP = length(pNormList);
    results_matrix_ADMM = zeros(num_pairs, 3, numP);

    %% call all methods

    for ind_pair = 1 : num_pairs
        frame1 = ind_pair;
        frame2 = frame1 + frame_gap;

        %% main loop
        pFs = [frame1, frame2]; % frame index
        parKnl = st('alg', 'cmum'); % type of affinity: only edge distance

        % algorithm parameter
        [pars, algs] = gmPar(2);

        % src
        wsSrc = cmumAsgSrc(tag, pFs, nIn, 'svL', 1);
        asgT = wsSrc.asgT;

        % feature
        parG = st('link', 'del'); % Delaunay triangulation for computing the graphs
        parF = st('smp', 'n', 'nBinT', 4, 'nBinR', 3); % not used, ignore it
        wsFeat = cmumAsgFeat(wsSrc, parG, parF, 'svL', 1);
        [gphs, XPs, Fs] = stFld(wsFeat, 'gphs', 'XPs', 'Fs');

        % affinity
        [KP, KQ] = conKnlGphPQU(gphs, parKnl);
        K = conKnlGphKU(KP, KQ, gphs);
        Ct = ones(size(KP));
        
        % undirected graph -> directed graph (for FGM-D)
        gphDs = gphU2Ds(gphs);
        KQD = [KQ, KQ; KQ, KQ];

        % setup the Laplacian and constraints
        [n1, n2] = size(Ct); 
        D_matrix = diag(sum(K)); d_vector = diag(D_matrix); % K is the M matrix in our manuscript, see Eq. (25)
        L_matrix = D_matrix - K; 
        A_1 = kron(ones(1, n2), sparse(eye(n1))); 
        A_2 = kron(sparse(eye(n2)), ones(1, n1));
        A = [A_1; A_2]; % one-to-one assignment constraint
        b = ones(size(A,1), 1);

        %% Apply All Methods
        % [1]. PM
        tic; asgPm = pm(K, KQ, gphs, asgT); asgPm.tim = toc; 
        results_matrix_PM(ind_pair, :) = [asgPm.acc, asgPm.obj , asgPm.tim];

        % [2]. SM
        asgSm = gm(K, Ct, asgT, pars{3}{:});
        results_matrix_SM(ind_pair, :) = [asgSm.acc, asgSm.obj , asgSm.tim];

        % [3]. IPFP-U
        asgIpfpU = gm(K, Ct, asgT, pars{5}{:}); 
        results_matrix_IPFP_U(ind_pair, :) = [asgIpfpU.acc, asgIpfpU.obj , asgIpfpU.tim];

        % [4]. IPFP-S
        asgIpfpS = gm(K, Ct, asgT, pars{6}{:});
        results_matrix_IPFP_S(ind_pair, :) = [asgIpfpS.acc, asgIpfpS.obj , asgIpfpS.tim];

        % [5]. RRWM
        asgRrwm = gm(K, Ct, asgT, pars{7}{:});
        results_matrix_RRWM(ind_pair, :) = [asgRrwm.acc, asgRrwm.obj , asgRrwm.tim];

        % [6]. FGM-U
        tic; asgFgmU = fgmU(KP, KQ, Ct, gphs, asgT, pars{8}{:}); asgFgmU.tim = toc; 
        results_matrix_FGM_U(ind_pair, :) = [asgFgmU.acc, asgFgmU.obj , asgFgmU.tim];

        % [7]. FGM-D
        tic; asgFgmD = fgmD(KP, KQD, Ct, gphDs, asgT, pars{9}{:}); asgFgmD.tim = toc; 
        results_matrix_FGM_D(ind_pair, :) = [asgFgmD.acc, asgFgmD.obj , asgFgmD.tim];


        % [8]. QP, quadprog

        % parameter setting 
        % n1 = nIn(1); n2 = nIn(2); 
        [n1, n2] = size(Ct); 
        D_matrix = diag(sum(K)); d_vector = diag(D_matrix); % K is the M matrix in our draft, see Eq. 67
        L_matrix = D_matrix - K; 
        A_1 = kron(ones(1, n2), eye(n1)); 
        A_2 = kron(eye(n2), ones(1, n1));
        A = [A_1; A_2]; % one-to-one assignment constraint
        b = ones(size(A,1), 1);
        lb = zeros(n1*n2, 1); ub = ones(n1*n2, 1); 
        % min_x   x' * L_matrix * x - d_vector' * x,  s.t. x \in [0,1], A x <= b. 

        % call quadprog
        % min_x  1/2 * x' * H * x + f' * x;   H = 2.*L_matrix, f = -d_vector; s.t. x \in [0,1],
        tic; [x_vec_QP, obj_QP_continuous, exitflag_QP,output_QP] = quadprog(2.*L_matrix, -d_vector, A, b, [], [],lb,ub); 

        % transform the continuous results to binary matching, through Hugarian algorithm
        X_QP = reshape(x_vec_QP, n1, n2);
        [C_QP,T_QP]=hungarian(1-X_QP);
        % x_vec_QP' * K * x_vec_QP
        % x_vec_QP' * L_matrix * x_vec_QP - d_vector' * x_vec_QP
        time_QP = toc;

        % evaluation
        acc_QP = matchAsg(C_QP, asgT);
        obj_QP_disc = vec(C_QP)' * K * vec(C_QP); 
        results_matrix_QP(ind_pair, :) = [acc_QP, obj_QP_disc, time_QP];
        
        
        % [9]. ADMM
        % min_x   x' * L_matrix * x - d_vector' * x,  s.t. x \in {0,1}, A x <= b.  

        for pID = 1:numP
                p = pNormList(pID); 

                acc_ADMM = 0;
                obj_ADMM_disc = 0;
                time_ADMM = 0;

                x0 = vec(asgSm.X);   % initialization by the result of SM, you can try the results of other methods
                params = struct('is_verbose',false,'std_threshold',1e-6,...
                    'stop_threshold',1e-6,'gamma_val',1.0,'gamma_factor',1,'initial_rho',20,...
                    'x0',x0,'learning_fact',1+3/100,'history_size',3,'rel_tol',1e-5,'max_iters',5e2, 'projection_lp', p);

                if num_node_remove~=0
                    % the constraint is an inequality    
                    [~,x_sol,obj_list,~,~,~,~,time_ADMM] = ADMM_bqp_linear_inequality(L_matrix,-d_vector,A,b,params);            
                else
                    % the constraint is an equality
                    [~,x_sol,obj_list,~,~,~,time_ADMM] = ADMM_bqp_linear_eq(L_matrix,-d_vector,A,b,params);            
                end
                figure; plot(1:length(obj_list), obj_list, '<-g')

                % evaluate the solution
                bin_x_sol = x_sol>=0.5;
                acc_ADMM = matchAsg(bin_x_sol, asgT);
                obj_ADMM_disc = vec(bin_x_sol)'*K*vec(bin_x_sol); 

                results_matrix_ADMM(ind_pair, :, pID) = [acc_ADMM, obj_ADMM_disc, time_ADMM];
        end
        
    end

    %% show the mean of results of all methods
    all_results = cat(3,results_matrix_PM,results_matrix_SM,results_matrix_IPFP_U,results_matrix_IPFP_S,...
        results_matrix_RRWM,results_matrix_FGM_U,results_matrix_FGM_D,results_matrix_QP,results_matrix_ADMM);

    m = squeeze(median(all_results,1));
    method_names = {'PM', 'SM', 'IPFP_U', 'IPFP_S', 'RRWM', 'FGM_U', 'FGM_D', 'QP', 'ADMM-L05', 'ADMM-L1', 'ADMM-L2', 'ADMM-L5', 'ADMM-L10', };
    
    
    all_results_struct = struct('all_results', all_results, 'params_ADMM', params); 
    all_results_struct.method_names = method_names;
    
    savePath = './matching/results/';
    if ~isdir(savePath), mkdir(savePath); end
    saveName = ['all_results_framegap_', num2str(frame_gap), '_nodeRemove_', num2str(num_node_remove), '_', generate_clock_str(), '.mat'  ];
    save([savePath, saveName], 'all_results_struct')

end
