%% Loading functions and Data generation
clear;
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Parameters for the data generating mechanism
rng(12242,'twister');

% load parameters for the model
run('./Parameters_setup_ground_truth.m')

num_dense=40; % Number of grids
Aval = 20;
A = diag([Aval, Aval, 750]); %<--- 

run('Parameters_setup_experiment.m')

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

%%
num_sim = 10;
%% Random design
design = 0; % 0: random design; 1: optimal design

for randomseed = 1:num_sim    
    rng(randomseed,'twister');

    flnm=strcat('../../Data/sim-results/A', num2str(Aval), 'Design', num2str(design),...
        'Seed',num2str(randomseed),'.mat'); 
    % Run analysis and design
    
    % Initialize starting values
    output= struct([]);
    for j = 1:num_threshold
        output(j).alpha = .1*ones(K_z+1,1);
        output(j).mu = zeros(K_z+1,1);
        output(j).s_sq = ones(K_z+1,1);
    end
    
    X_g = zeros((N-1)*num_trials_batch + num_trials_first,K_z);
    locations_trials = zeros( (N-1)*num_trials_batch + num_trials_first,num_sources);
    Y_g = zeros((N-1)*num_trials_batch + num_trials_first,num_threshold);
    
    % Initialize storage:
    output_random =cell(N/10,1);
    time_record = zeros(N,1);
    
    % Conducting experiments
    for outer_i = 1:(N/10)
        t_start = (outer_i-1)*10+1;
        t_end = outer_i*10;
        run('Experiment_full.m');
        output_random{outer_i}=output;
    end
    X_random = X_g;
    mpp_random = mpp_new;
    time_record_random = time_record;
    location_random = locations_trials;
    
    % Evaluate the performance:
    % 1, normalized reconstruction error of connectivity
    % 2, AUC for connectivity
    % 3, normalized amplitudes reconstruction
    output_eva = output_random;
    run('Experiment_evaluate.m')
    
    NRE_conn_random = NRE_conn;
    NRE_amp_random = NRE_amp;
    AUC_conn_random = AUC_conn;

    save(flnm, 'output_random','X_random',...
        'location_random','mpp_random','time_record_random',...
        'NRE_conn_random','NRE_amp_random','AUC_conn_random')
end




%% Optimal design

design = 1;

for randomseed = 1:num_sim
    rng(randomseed,'twister');

    flnm=strcat('../../Data/sim-results/A', num2str(Aval), 'Design', num2str(design),...
        'Seed',num2str(randomseed),'.mat');
    % Initialize starting values
    output= struct([]);
    for j = 1:num_threshold
        output(j).alpha = .1*ones(K_z+1,1);
        output(j).mu = zeros(K_z+1,1);
        output(j).s_sq = ones(K_z+1,1);
    end
    
    
    X_g = zeros((N-1)*num_trials_batch + num_trials_first,K_z);
    locations_trials = zeros( (N-1)*num_trials_batch + num_trials_first,num_sources);
    Y_g = zeros((N-1)*num_trials_batch + num_trials_first,num_threshold);
    
    % Initialize storage:
    output_optimal =cell(N/10,1);
    time_record = zeros(N,1);
    
    % Conducting experiments
    for outer_i = 1:(N/10)
        t_start = (outer_i-1)*10+1;
        t_end = outer_i*10;
        run('Experiment_full.m');
        output_optimal{outer_i}=output;
    end
    X_optimal= X_g;
    mpp_optimal = mpp_new;
    time_record_optimal = time_record;
    location_optimal = locations_trials;
    
    % Evaluate the performance:
    % 1, normalized reconstruction error of connectivity
    % 2, AUC for connectivity
    % 3, normalized amplitudes reconstruction
    output_eva = output_optimal;
    run('Experiment_evaluate.m')
    
    NRE_conn_optimal = NRE_conn;
    NRE_amp_optimal = NRE_amp;
    AUC_conn_optimal = AUC_conn;
    
    save(flnm, 'output_optimal','X_optimal',...
        'location_optimal','mpp_optimal','time_record_optimal',...
        'NRE_conn_optimal','NRE_amp_optimal','AUC_conn_optimal')
    
end
