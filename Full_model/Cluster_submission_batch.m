% Loading functions and Data generation
tic
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Paramters in the simulations
num_dense_grid          = [80   80  80  80  80  80  80  80  80  80  80  80  80];
V_homo_grid             = [0    0   0   0   0   0   0   0   0   0   0   0   0];

Aval_grid               = [400  100 900 400 400 400 400 400 400 400 400 400 400]; 
freq_pen_grid           = [1    1   1   1.5 2   1   1   1   1   1   1   1   1];
num_trials_batch_grid   = [20   20  20  20  20  40  100 20  20  20  20  20  20]; 
layer_feature_grid      = [0    0   0   0   0   0   0   1   0   0   0   0   0];
trials_specific_var_grid= [0    0   0   0   0   0   0   0   1   0   0   0   0];
random_prop_grid        = [0    0   0   0   0   0   0   0   0   0.2 0.4 0   0.2];
num_random_grid         = [4    4   4   4   4   4   4   4   4   4   4   2   2]; % number of initial trials

%%
for i_setting = 1:length(num_dense_grid)
    
    %---------------------------------------------------------------------%
    num_dense= num_dense_grid(i_setting); % Number of grids
    Aval = Aval_grid(i_setting); 
    V_homo=V_homo_grid(i_setting);
    freq_pen = freq_pen_grid(i_setting);
    num_random =num_random_grid(i_setting);
    trials_specific_variance= trials_specific_var_grid(i_setting);
    random_prop = random_prop_grid(i_setting);
    num_trials_batch=num_trials_batch_grid(i_setting);
    layer_feature=layer_feature_grid(i_setting);
    
    length_memory = num_trials_batch;
    num_peaks = 2*num_trials_batch;
    gamma_threshold = 0.1; % For the sparse EM algorithm 
    
    num_samples=50; % Number of samples to estimate the expected entropies
    num_sources = 4;  % Number of locations to stimulate in each trial
    N=4000/num_trials_batch; % Number of batches
    
    
    cell_feature_priors=struct();
    num_layers = 7;
    if V_homo == 1
        cell_feature_priors.Vthre_mean = [0  5  5  5  5  5 5]; % gaussian
        cell_feature_priors.Vreset_mean = [0  -50  -50  -50  -50  -50  -50]; % gaussian
        
    else
        cell_feature_priors.Vthre_mean = [0  3  3  5  5  7 7]; % gaussian
        cell_feature_priors.Vreset_mean = [0  -45  -45  -50  -50  -55  -55]; % gaussian
    end
    cell_feature_priors.Vthre_std = 1* ones(num_layers,1); % gaussian
    cell_feature_priors.Vreset_std =  1* ones(num_layers,1); % gaussian
    
    %---------------------------------------------------------------------%
    
    % Parameters for the data generating mechanism
    rng(12242,'twister');
    % load parameters for the model
    
    A = diag([Aval, Aval, 1500]); 
    run('./Data_generation/Parameters_setup_LIF.m')
    num_I_Stim=1; % set the stimulus to be a constant value, e.g., 100
    I_e_vect=[0;I_e(:,num_I_Stim)];
    evoked_params.stim_start = min(find(I_e_vect>10));
    evoked_params.stim_end = max(find(I_e_vect>10));
    
    % Stimulation:
    num_sources = 4;  % Number of locations to stimulate in each trial
    grid_index = 1:size(pi_dense_all,2);
    % n_cell = size(pi_dense_all,1);
    % n_cell_local = size(pi_dense_local,1);
    
    k_minimum = 0.001; % minimum stimulus intensity to consider
    k_offset = 0.04; % the spatial mark = k_offset*spatial distance etc
    
    % Stimulation:
    % Set seed:
    rng(arg1,'twister');
    
    % Parameters
    %
    sqrt_transform = false; % whether to use squared transformation
    
    % Parameters for the working model
    num_threshold=10; % number of bins to use
    mark = 0; % 0: amplitude; 1: latency.
    obj_function = @joint_sparsity_weight_entropy; %
    
    
    t_vect=0:data_params.dt:data_params.T;
    V_vect=zeros(1,length(t_vect));
    num_trials_first =max(200, ceil(n_cell_local/num_sources)); % Number of trials in the first batches
    
    %---------------------------------------------------------------------%
    % Designing stage
    
    % Design
    %design = 0
    % Run analysis and design
    % Initialize starting values
    output= struct([]);
    for j = 1:num_threshold
        output(j).alpha = .1*ones(n_cell_local+1,1);
        output(j).mu = zeros(n_cell_local+1,1);
        output(j).s_sq = ones(n_cell_local+1,1);
        output(j).threshold = [];
    end
    output_ini = output;
    X_g = zeros((N-1)*num_trials_batch + num_trials_first,n_cell_local);
    locations_trials = zeros( (N-1)*num_trials_batch + num_trials_first,num_sources);
    Y_g = zeros((N-1)*num_trials_batch + num_trials_first,num_threshold);
    % Initialize storage:
    output_crude =cell(N,1);
    time_record = zeros(N,1);
    
    % Conducting experiments
    for t = 1:N
        run('./Design/Design_LIF.m');
        output(num_threshold).threshold= threshold;
        output_crude{t}=output;
    end
    stimuli_size_local = k_offset.*X_g;
    
   
    
    
    
    %------------------------------------%
    % Estimating the marginal firing rate 
    sigma_unknown=1;
    I_stimuli = I_e_vect;
    T=75;
    dt=1;
    t_vect=0:dt:T;
    if layer_feature == 1
        V_thresholds = local_V_th_layer;
        V_resets = local_V_reset_layer;
    else
        V_thresholds = local_V_th;
        V_resets = local_V_reset;
    end
    
    n_stimuli_grid=10;
    n_trial = size(stimuli_size_local,1);
    n_cell_local = size(stimuli_size_local,2);
    % Prepare the batch data sets
    estimated_intensity=cell(n_trial, n_cell_local);
    evoked_cell = cell(n_trial,1);
    
    stimuli_temp = stimuli_size_local;
    run('./Inference/Expected_intensity_v3.m');
    estimated_intensity=M_intensity;
    
    for i_trial = 1:n_trial
        evoked_cell_index = 0; % 0: background evnets
        for i_cell = 1:n_cell_local
            k = stimuli_size_local(i_trial, i_cell);
            if k > k_minimum
                evoked_cell_index = [evoked_cell_index i_cell];
            end
        end
        evoked_cell{i_trial} = evoked_cell_index;
    end
    expected_all = zeros(n_trial,n_cell_local);
    for i_trial = 1:n_trial
        for i_cell = 1:n_cell_local
            expected_all(i_trial,i_cell)=sum(estimated_intensity{i_trial, i_cell} );
        end
    end
    %---------------------------------------------------------------------%
    
    
    convergence_epsilon = 0.01;
    maxit = 100;
    n_gibbs_sample = 100;
    n_burnin = 100;
    n_skip = 5;
    f_background = bg_params.firing_rate/1000;
    mean_background = bg_params.mean;
    sigma_background  = bg_params.sigma;
    
    
    t_seq = [1 2 3 4 5*(1: N/5)];
    
    output_EM = cell(length(t_seq),1);
    output_Gibbs = cell(length(t_seq),1);
    output_sparse = cell(length(t_seq),1);
    for ind_t = 1:length(t_seq)
        t=t_seq(ind_t); % iterate:
        fprintf('Batch %d\n', t);
        % New data in this batch
        
        idx_trials_this_batch =  1: (num_trials_first + ...
            (t - 1)*num_trials_batch );
        n_trial_update = max(idx_trials_this_batch);
        n_trial_temp = max(idx_trials_this_batch);
        
        
        
        %---------------------------------------------------------------------%
    % Initialization of the full model 
        output= output_crude{t};
        overall_connectivity = zeros(n_cell_local,1);
        overall_mark = zeros(n_cell_local,1);
        for j = 1:num_threshold
            overall_connectivity = max(overall_connectivity, ...
                output(j).alpha(2:end));
            for i_cell = 1:n_cell_local
                if overall_connectivity(i_cell) == output(j).alpha(i_cell+1)
                    overall_mark(i_cell)= (output(num_threshold).threshold(j)+...
                        output(num_threshold).threshold(j+1))/2;
                end
            end
        end
        
        
        %----------------------------------------------%
        % EM
        gamma_old = overall_connectivity;
        mu_old = overall_mark;
        sigma_old = 1*ones(n_cell_local,1);
        
        sparsity =0;
        tstart = toc;
        run('./Inference/Simulation_EM_batch.m');
        tend=toc;
        output_EM{ind_t}.sigma_samples = sigma_samples(:,end);
        output_EM{ind_t}.gamma_samples = gamma_samples(:,end);
        output_EM{ind_t}.mu_samples = mu_samples(:,end);
        %output_EM{ind_t}.soft_assignments_samples= soft_assignments_samples;
        output_EM{ind_t}.delta_t = tend-tstart;
        %----------------------------------------------%
        
        %----------------------------------------------%
        % EM with sparsity 
        gamma_old = overall_connectivity;
        mu_old = overall_mark;
        sigma_old = 1*ones(n_cell_local,1);
        
        sparsity =1;
        tstart = toc;
        run('./Inference/Simulation_EM_batch.m');
        tend=toc;
        output_sparse{ind_t}.sigma_samples = sigma_samples(:,end);
        output_sparse{ind_t}.gamma_samples = gamma_samples(:,end);
        output_sparse{ind_t}.mu_samples = mu_samples(:,end);
        output_sparse{ind_t}.delta_t = tend-tstart;
        %
        
        %--------------------------------------------------------%
        % Gibbs
        gamma_current = overall_connectivity;
        mu_current = overall_mark;
        sigma_current =4*ones(n_cell_local,1);
        %
        mu_m_hyper =overall_mark;
        mu_v_hyper = ones(n_cell_local,1); % inverse of prior variance
        sigma_alpha = ones(n_cell_local,1);
        sigma_beta = ones(n_cell_local,1);
        
        gamma_alpha = ones(n_cell_local,1);
        gamma_beta = ones(n_cell_local,1);
        
        
        tstart = toc;
        run('./Inference/Simulation_Gibbs_batch.m');
        tend=toc;
        output_Gibbs{ind_t}.sigma_samples = sigma_samples;
        output_Gibbs{ind_t}.gamma_samples = gamma_samples;
        output_Gibbs{ind_t}.mu_samples = mu_samples;
        output_Gibbs{ind_t}.delta_t = tend-tstart;
    end
    
    % record data
    %
    flnm = strcat('./Data/Setting',num2str(i_setting),'Design', num2str(design));
    
    save(strcat(flnm, 'Rep',num2str(arg1),'_Gibbs.mat'),'output_Gibbs');
    save(strcat(flnm, 'Rep',num2str(arg1),'_Sparse.mat'),'output_Sparse');
    save(strcat(flnm, 'Rep',num2str(arg1),'_EM.mat'),'output_EM');
    save(strcat(flnm, 'Rep',num2str(arg1),'_Data.mat'),'mpp_new','stimuli_size_local');
    save(strcat(flnm, 'Rep',num2str(arg1),'_Crude.mat'),'output_crude','time_record','threshold');
    save(strcat(flnm, 'Rep',num2str(arg1),'_Trials.mat'),'locations_trials');
    
    save(strcat(flnm, 'Rep',num2str(arg1),'_Truth.mat'),...
    'bg_params','neuron_features','neuron_locations','num_threshold', ...
        'time_record', 'k_minimum', 'N','num_trials_first','num_trials_batch',...
        'local_locations','local_amplitudes','evoked_params', 'local_V_th','local_V_reset',...
        'local_V_th_layer','local_V_reset_layer',...
        'local_sigma_within','local_sigma_across', 'local_gamma','Z');
    
    
end