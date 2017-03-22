%% Loading functions and Data generation
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Parameters for the data generating mechanism
rng(12242,'twister');
% load parameters for the model
% run('./Parameters_setup_ground_truth.m')
run('./Data_generation/Parameters_setup_LIF.m')

num_dense=100; % Number of grids
Aval = 400; % Very high precision 64, standard 400;
A = diag([Aval, Aval, 1500]); %<---
run('./Data_generation/Parameters_setup_experiment.m')
num_I_Stim=1; % set the stimulus to be a constant value, e.g., 100
%
k_basic = 0.04;
n_cell = length(all_amplitudes);
num_sources = 4;  % Number of locations to stimulate in each trial
%% Stimulation:
grid_index = 1:size(pi_dense_all,2);
n_cell_local = size(pi_dense_local,1);
k_minimum = 0.001; % minimum stimulus intensity to consider

%% Parameters
%
N=200; % Number of batches
num_trials_first = 200; % Number of trials in the first batches

num_samples=50; % Number of samples to estimate the expected entropies
num_trials_batch=20;
num_sources = 4;  % Number of locations to stimulate in each trial
num_peaks = 20;


sqrt_transform = false; % whether to use squared transformation
num_threshold=10; % number of bins to use
mark = 0; % 0: amplitude; 1: latency.
obj_function = @joint_sparsity_weight_entropy; %

% Fix the stimulation to be 100 for now
num_I_Stim=1; % set the stimulus to be a constant value, e.g., 100
% 1: 100; 6: 50; 11: 25;

%% Random design
if design ==0
        % Run analysis and design
        
        % Initialize starting values
        output= struct([]);
        for j = 1:num_threshold
            output(j).alpha = .1*ones(n_cell_local+1,1);
            output(j).mu = zeros(n_cell_local+1,1);
            output(j).s_sq = ones(n_cell_local+1,1);
        end
        
        X_g = zeros((N-1)*num_trials_batch + num_trials_first,n_cell_local);
        locations_trials = zeros( (N-1)*num_trials_batch + num_trials_first,num_sources);
        Y_g = zeros((N-1)*num_trials_batch + num_trials_first,num_threshold);
        
        % Initialize storage:
        output_random =cell(N/10,1);
        time_record = zeros(N,1);
        
        % Conducting experiments
        for outer_i = 1:(N/10)
            t_start = (outer_i-1)*10+1;
            t_end = outer_i*10;
            %run('Experiment_full.m');
            run('./Design/Design_LIF.m');
            output_random{outer_i}=output;
        end
        
        stimuli_size_local = k_basic.*X_g;
        % mpp_new;
        % time_record;
        %locations_trials;
        
        % Evaluate the performance:
        % 1, normalized reconstruction error of connectivity
        % 2, AUC for connectivity
        % 3, normalized amplitudes reconstruction
        %output_eva = output_random;
        %run('Experiment_evaluate.m')
        
        %NRE_conn_random = NRE_conn;
        %NRE_mark_random = NRE_mark;
        %AUC_conn_random = AUC_conn;
        
        
       flnm = strcat('./Data/RandomDesign.mat');
        save(flnm,'bg_params','neuron_features','neuron_locations', ...
            'time_record', 'k_minimum',...
        'local_locations','local_amplitudes','evoked_params', 'local_V_th','local_V_reset','local_E_L', ...
            'Z');
    
       flnm = strcat('./Data/RandomDesign_data.mat');
        save(flnm,'mpp_new','stimuli_size_local');
        
        
elseif design == 1 % Optimal design
        
        % Initialize starting values
        output= struct([]);
        for j = 1:num_threshold
            output(j).alpha = .1*ones(n_cell_local+1,1);
            output(j).mu = zeros(n_cell_local+1,1);
            output(j).s_sq = ones(n_cell_local+1,1);
        end
        X_g = zeros((N-1)*num_trials_batch + num_trials_first,n_cell_local);
        locations_trials = zeros( (N-1)*num_trials_batch + num_trials_first,num_sources);
        Y_g = zeros((N-1)*num_trials_batch + num_trials_first,num_threshold);
       
        % Initialize storage:
        output_optimal =cell(N/10,1);
        time_record = zeros(N,1);
        
        % Conducting experiments
        for outer_i = 1:(N/10)
            t_start = (outer_i-1)*10+1;
            t_end = outer_i*10;
            run('./Design/Design_LIF.m');
            output_optimal{outer_i}=output;
        end
        stimuli_size_local = k_basic.*X_g;
        % mpp_new;
        %time_record;
        % locations_trials;
        
        % Evaluate the performance:
        % 1, normalized reconstruction error of connectivity
        % 2, AUC for connectivity
        % 3, normalized amplitudes reconstruction
%         output_eva = output_optimal;
%         run('Experiment_evaluate.m')
%         
%         NRE_conn_optimal = NRE_conn;
%         NRE_mark_optimal = NRE_mark;
%         AUC_conn_optimal = AUC_conn;
%         
        
        
       flnm = strcat('./Data/OptimalDesign.mat');
        save(flnm,'bg_params','neuron_features','neuron_locations', ...
           'time_record', 'k_minimum',...
        'local_locations','local_amplitudes','evoked_params', 'local_V_th','local_V_reset','local_E_L', ...
            'Z');
    
       flnm = strcat('./Data/OptimalDesign_data.mat');
        save(flnm,'mpp_new','stimuli_size_local');
        
    
end