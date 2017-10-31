% cell_parameters_type=1;
% prior_info_type=1;
% model_type=1;
% design_type =1;
% assignment_type=1;


%% Generate cellular parameters
% tuning parameter: cell_parameters_type;
run('./Simulation_parameters.m')
%% Preprocessing
run('./Simulation_preprocessing.m')
%% Designing experiment
%-------------------------------------------------%
%% Parameters in the design stage
n_replicates=1; % number of replicates for each trial
n_spots_per_trial = 4;trial_max=2000;
K_connected=10;K_undefined=20;

switch assignment_type
    case 1
        disconnected_threshold = 0.1;connected_threshold = 0.4;
        quantiles_prob=[0.1;0.9];
    case 2
        disconnected_threshold = 0.05;connected_threshold = 0.95;
        gamma_bound.up=1;gamma_bound.low=0.05;
        quantiles_prob=[0.1;0.90];
    otherwise
end
        
% if the postsynaptic cell has good connection, switch to 
disconnected_threshold_good = 0.2;connected_threshold_good = 0.5;
quantiles_prob_good=[0.1;0.90];

multiplier_connected=1;
single_spot_threshold=12; %12*10/4=30 (only 30 trials on these cells)
% switch to single spot stimulation (this can be a function of n_spots_per_trial
multispot_change_threshold=0.05; % switch if there aren't enough changes 

% Define tuning parameters in the VI
maxit=1000;
S=50;epsilon=0.01;eta=1;eta_max=2;
background_rt=background_rate*time_max; % raw probability of firing within a trial
n_MC_samples=50;


num_trials_per_batch=200;
% tuning parameter: prior_info_type;
run('./Simulation_tuning.m')

%%

disconnected_params=struct([]);
disconnected_params(1).prop=0.8;
disconnected_params(1).alpha=-4;
disconnected_params(1).beta=0;

regular_params=struct([]);
regular_params(1).prop=0.2;
regular_params(1).alpha=0.6;
regular_params(1).beta=0;

som_params=struct([]);
som_params(1).prop=0.2;
som_params(1).alpha=-1.4;
som_params(1).beta=0;




% tuning parameters: model_type, design_type, assignment_type
run('./Simulation_experiment.m')

