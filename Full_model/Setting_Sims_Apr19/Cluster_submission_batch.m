% Loading functions and Data generation
tic
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Paramters in the simulations
num_dense_grid = [80 40 120 80 80 80 80 80 80 80]; % desnity of grids
Aval_grid = [400 400 400 100 900 400 400 400 400 400]; % precision of stimuli
num_random_grid = [4 4 4 4 4 2 7 5 5 4]; % number of initial trials
V_homo_grid = [1*ones(1,9) 0];
freq_pen_grid = [1*ones(1,7) 1.2 1.5 0];
%%
for i_setting = 1:length(num_dense_grid)
    num_dense= num_dense_grid(i_setting); % Number of grids
    Aval = Aval_grid(i_setting); % Very high precision 64, standard 400;
    num_random =num_random_grid(i_setting);
    V_homo=V_homo_grid(i_setting);
    freq_pen = freq_pen_grid(i_setting);
    
    length_memory = 20;
    
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
    cell_feature_priors.Vreset_std =  2* ones(num_layers,1); % gaussian
    
    
    
    % Parameters for the data generating mechanism
    rng(12242,'twister');
    % load parameters for the model
    
    A = diag([Aval, Aval, 1500]); %<---
    
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
    N=200; % Number of batches
    
    num_samples=50; % Number of samples to estimate the expected entropies
    num_trials_batch=20;
    num_sources = 4;  % Number of locations to stimulate in each trial
    num_peaks = 40;
    num_trials_first =max(200, ceil(n_cell_local/num_sources)); % Number of trials in the first batches
    sqrt_transform = false; % whether to use squared transformation
    
    % Parameters for the working model
    num_threshold=10; % number of bins to use
    mark = 0; % 0: amplitude; 1: latency.
    obj_function = @joint_sparsity_weight_entropy; %
    
    
    t_vect=0:data_params.dt:data_params.T;
    V_vect=zeros(1,length(t_vect));
    
    
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
    
   
    
    
    %
    %------------------------------------%
    % Inference
    
    convergence_epsilon = 0.01;
    maxit = 100;
    
    sigma_unknown=1;
    
    I_stimuli = I_e_vect;
    
    T=75;
    dt=1;
    t_vect=0:dt:T;
    
    V_thresholds = local_V_th;
    V_resets = local_V_reset;
    
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
    
    
    t_seq = [1 2 3 4 5 5*(1: N/5)];
    
    output_EM = cell(length(t_seq),1);
    output_working = cell(length(t_seq),1);
    for ind_t = 1:length(t_seq)
        t=t_seq(ind_t); % iterate:
        fprintf('Batch %d\n', t);
        % New data in this batch
        
        idx_trials_this_batch =  1: (num_trials_first + ...
            (t - 1)*num_trials_batch );
        
        n_trial_update = max(idx_trials_this_batch);
        n_trial_temp = max(idx_trials_this_batch);
        
        % Initialize the estimates
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
        
        gamma_old = overall_connectivity;
        mu_old = overall_mark;
        sigma_old = 1*ones(n_cell_local,1);
        
        
        f_background = bg_params.firing_rate/1000;
        mean_background = bg_params.mean;
        sigma_background  = bg_params.sigma;
        
        
        tstart = toc;
        run('./Inference/Simulation_EM_batch.m');
        tend=toc;
        output_EM{ind_t}.sigma_samples = sigma_samples(:,end);
        output_EM{ind_t}.gamma_samples = gamma_samples(:,end);
        output_EM{ind_t}.mu_samples = mu_samples(:,end);
        %output_EM{ind_t}.soft_assignments_samples= soft_assignments_samples;
        output_EM{ind_t}.delta_t = tend-tstart;
        %
        
        %     gamma_current = overall_connectivity;
        %     mu_current = overall_mark;
        %     sigma_current =4*ones(n_cell_local,1);
        %     %
        %     mu_m_hyper =overall_mark;
        %     mu_v_hyper = ones(n_cell_local,1); % inverse of prior variance
        %     sigma_alpha = ones(n_cell_local,1);
        %     sigma_beta = ones(n_cell_local,1);
        %
        %     gamma_alpha = ones(n_cell_local,1);
        %     gamma_beta = ones(n_cell_local,1);
        %
        %     f_background = bg_params.firing_rate/1000;
        %     mean_background = bg_params.mean;
        %     sigma_background  = bg_params.sigma;
        %
        %     tstart = toc;
        %     run('./Inference/Simulation_integral_batch.m');
        %     tend=toc;
        %     output_Gibbs{t}.sigma_samples = sigma_samples;
        %     output_Gibbs{t}.gamma_samples = gamma_samples;
        %     output_Gibbs{t}.mu_samples = mu_samples;
        %     output_Gibbs{t}.soft_assignments_samples= soft_assignments_samples;
        %     output_Gibbs{t}.delta_t = tend-tstart;
        %
        %
    end
    
    % record data
    %
    flnm = strcat('./Data/Design', num2str(design),'Setting',num2str(i_setting));
    
    save(strcat(flnm, 'Rep',num2str(arg1),'_Working.mat'),'output_working');
    
    save(strcat(flnm, 'Rep',num2str(arg1),'_EM.mat'),'output_EM');
    
    save(strcat(flnm, 'Rep',num2str(arg1),'_Truth.mat'),...
    'bg_params','neuron_features','neuron_locations','num_threshold', ...
        'time_record', 'k_minimum', 'N','num_trials_first','num_trials_batch',...
        'local_locations','local_amplitudes','evoked_params', 'local_V_th','local_V_reset',...
        'local_sigma', 'local_gamma','Z');
    
    % flnm=strcat('./Data/Design',num2str(design),'Batch_Gibbs_Set',num2str(arg1),'.mat');
    % save(flnm,'output_Gibbs');
    %flnm=strcat('./Data/Design',num2str(design),'Batch_Intensity_Set',num2str(arg1),'.mat');
    %save(flnm,'estimated_intensity');
   
    save(strcat(flnm, 'Rep',num2str(arg1),'_Data.mat'),'mpp_new','stimuli_size_local');
    save(strcat(flnm, 'Rep',num2str(arg1),'_Crude.mat'),'output_crude','time_record','threshold');
    save(strcat(flnm, 'Rep',num2str(arg1),'_Trials.mat'),'locations_trials');
    
end

%%

%      figure(ind_t)
%      plot(local_gamma,output_EM{ind_t}.gamma_samples,'.','col',[0,1,0,1], 'markersize', 10);
%      hold on;
%      %plot(jittered_gamma,output_working{ind_t}.mean_gamma,'.','col',[0,0,1,1], 'markersize', 10);
%      %hold on;
%
%     % plot(jittered_gamma,output_crude{t_seq(ind_t)}.mean_gamma,'.','col',[1,0,0,1], 'markersize', 10);
%      %hold on;
%
%      line([0 1],[0 1]);
%
%      xlim([-0.1,1.1]);
%      ylim([-0.1,1.1]);
%
%
%      %%
%      plot(local_amplitudes,output_EM{ind_t}.mu_samples,'.','col',[0,1,0,1], 'markersize', 10);
%      hold on;
%         line([0 12],[0 12]);
%