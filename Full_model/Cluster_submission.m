% Loading functions and Data generation
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));

% Generate neurons:
rng(12242,'twister');
run('./Data_generation/Parameters_setup_LIF.m')
num_dense=100; % Number of grids
Aval = 400; % Very high precision 64, standard 400;
A = diag([Aval, Aval, 1500]); %<---
run('./Data_generation/Parameters_setup_experiment.m')
num_I_Stim=1; % set the stimulus to be a constant value, e.g., 100
%
n_trial = 4000;
k_basic = 0.04;
n_cell = length(all_amplitudes);
num_sources = 4;  % Number of locations to stimulate in each trial
%------------------------------------------------------------------%

% Stimulation:
% Set seed:
rng(arg1,'twister');

grid_index = 1:size(pi_dense_all,2);
K_z = size(pi_dense_local,1);
X_next = zeros(size(pi_dense_local,1), n_trial);
X_next_all = zeros(size(pi_dense_all,1), n_trial);
for l = 1:n_trial
    locations_next = randsample(grid_index, num_sources);
    locations_trials(l,:)=locations_next;
    X_next(:,l) = sum(pi_dense_local(:, locations_next) ,2);
    X_next_all(:,l) = sum(pi_dense_all(:,locations_next),2);
end
stimuli_size = k_basic*X_next_all';
stimuli_size_local = k_basic*X_next';

n_cell_local = size(stimuli_size_local,2);
k_minimum = 0.001; % minimum stimulus intensity to consider

% Generate data using the LIF-GLM model
% Note: 1 means using the LIF model.
exact_crossing = 0;
run('./Data_generation/Experiment_LIF.m');

flnm=strcat('./Data/example', num2str(arg1),'.mat');
save(flnm,'mpp_new','stimuli_size','stimuli_size_local');
flnm = strcat('./Data/truth', num2str(arg1),'.mat');
save(flnm,'bg_params','stimuli_size','stimuli_size_local','neuron_features','neuron_locations', 'k_minimum',...
    'local_locations','local_amplitudes','evoked_params', 'local_V_th','local_V_reset','local_E_L', ...
    'Z');


%------------------------------------%
% Inference
% Inference with the working model
num_threshold=10; % number of bins to use
mark = 0; % 0: amplitude; 1: latency.
n_trial = size(stimuli_size,1);
n_cell_local = size(stimuli_size_local,2);
tic
tstart = toc;
run('./Inference/Simulation_crude.m');
tend=toc;
t_delta = tend-tstart;

%------------------------%
flnm=strcat('./Data/Full_crude', num2str(arg1),'.mat');
save(flnm,'t_delta','overall_connectivity','overall_mark');
%------------------------%


% Run the full model
% 1. Estimate the `marginal' firing rate using paramters from LIF-GLM model
% Note: we use the true values in this simulation
load('./Environments/current_template.mat'); %Contains the vector norm_average_current
I_eg_100mV=abs(norm_average_current(1:20:20*75).*100);
I_eg_50mV=abs(norm_average_current(1:20:20*75).*50);
I_eg_25mV=abs(norm_average_current(1:20:20*75).*25);
I_e=[repmat(I_eg_100mV',1,5) repmat(I_eg_50mV',1,5) repmat(I_eg_25mV',1,5)];
num_I_Stim=1;
I_e_vect=[0;I_e(:,num_I_Stim)];
I_stimuli = I_e_vect;

% Stochastic components of voltages

T=75;
dt=1;
t_vect=0:dt:T;

V_thresholds = local_V_th;
V_resets = local_V_reset;
E_Ls = local_E_L; % need to know this a priori 

n_stimuli_grid=40;
k_basic = 0.04;
exact_crossing = 0;
run('./Inference/Expected_intensity_v2.m');

%--------------------%
flnm=strcat('./Data/estimated_intensity', num2str(arg1),'.mat');
save(flnm,'M_intensity');
%---------------------%

evoked_cell = cell(n_trial,1);
for i_trial = 1:n_trial
    evoked_cell_index = 0; % 0: background evnets
    for i_cell = 1:n_cell_local
        k = stimuli_size_local(i_trial, i_cell);
        if k > k_basic/10
            evoked_cell_index = [evoked_cell_index i_cell];
        end
    end
    evoked_cell{i_trial} = evoked_cell_index;
end

%
convergence_epsilon = 0.01;
maxit = 1000;

tic
n_trial_update = 4000;
sigma_unknown=1;
tstart=toc;
run('./Inference/Simulation_EM.m');
tend=toc;
t_delta = tend-tstart;

flnm=strcat('./Data/Full_EM', num2str(arg1),'.mat');
save(flnm,'t_delta','sigma_samples','gamma_samples','mu_samples', ...
    'soft_assignments_samples');

%
n_gibbs_sample = 400;
n_burnin = 800;
n_skip = 20;

tic
n_trial_update = 1000;
sigma_unknown=1;
tstart=toc;
run('./Inference/Simulation_integral.m');
tend=toc;
t_delta = tend-tstart;

flnm=strcat('./Data/Full_minibatch_int', num2str(arg1),'.mat');
save(flnm,'t_delta','sigma_samples','gamma_samples','mu_samples', ...
    'soft_assignments_samples');


