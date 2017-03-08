%% Loading functions and Data generation
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Parameters for the data generating mechanism
rng(12242,'twister');
% load parameters for the model
% run('./Parameters_setup_ground_truth.m')
run('./Data_generation/Parameters_setup_LIF.m')

num_dense=40; % Number of grids
Aval = 400;
A = diag([Aval, Aval, 1500]); %<--- 
run('./Data_generation/Parameters_setup_experiment.m')
num_I_Stim=1; % set the stimulus to be a constant value, e.g., 100
% 
n_trial = 4000;
k_basic = 0.1;
n_cell = length(all_amplitudes);
num_sources = 4;  % Number of locations to stimulate in each trial
%% Stimulation:
grid_index = 1:size(pi_dense_all,2);
K_z = size(pi_dense_local,1);

X_next = zeros(size(pi_dense_local,1), n_trial);
X_next_all = zeros(size(pi_dense_all,1), n_trial);
for l = 1:n_trial
    locations_next = randsample(grid_index, num_sources);
    locations_trials(l,:)=locations_next;
    X_next(:,l) = min(.95,sum(pi_dense_local(:, locations_next) ,2));
    X_next_all(:,l) = min(.95,sum(pi_dense_all(:,locations_next),2));
end

stimuli_size = k_basic*X_next_all';
stimuli_size_local = k_basic*X_next';


n_cell_local = size(stimuli_size_local,2);

flnm = strcat('./Data/truth.mat');
save(flnm,'stimuli_size','stimuli_size_local','neuron_features','neuron_locations',...
'local_locations','local_amplitudes','local_V_th','local_V_reset','local_E_L');
%% Generate data using the LIF-GLM model 
% Note: 1 means using the LIF model.
exact_crossing = 0;
%%
 run('./Data_generation/Experiment_LIF.m');
 
flnm=strcat('./Data/example1.mat');
save(flnm,'mpp_new','stimuli_size','stimuli_size_local');
