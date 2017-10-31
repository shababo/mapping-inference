addpath(genpath('../../../mapping-inference'));
%% Specify the setting in this simulation

% for i_sim = [12 15]
%     for i_seed = 1:20

cell_parameters_type=1;
prior_info_type=1;
model_type=1;
design_type =1;
assignment_type=1;
%%
switch cell_parameters_type
    case 1 % normal
        disp('normal gains and gammas')
    case 2 %
        disp('log-normal gains, normal gammas')
    case 3 %
        disp('normal gains, many weak gammas')
    otherwise
        disp('No inputs, use normal gains')
end

switch prior_info_type
    case 1
        disp('Use good prior')
    case 2 %
        disp('Use uninformative prior')
    case 3 %
        disp('Biased prior')
    otherwise
end

%% Generate cellular parameters
rng(i_seed,'twister');
% tuning parameter: cell_parameters_type;
run('./Simulation_parameters.m')
%% Preprocessing
run('./Simulation_preprocessing.m')
%% Designing experiment
%-------------------------------------------------%
%% Parameters in the design stage
n_replicates=1; % number of replicates for each trial
n_spots_per_trial = 4;trial_max=2000;
K_connected=10;
K_undefined=20;
assignment_type=1; % 1 gamma 2 spike 3 gain 4 quantiles of gamma

switch assignment_type
    case 1
        disconnected_threshold = 0.1;connected_threshold = 0.4;
        quantiles_prob=[0.05;0.95];
    case 2
        disconnected_threshold = 0.05;connected_threshold = 0.95;
        gamma_bound.up=1;gamma_bound.low=0.05;
        quantiles_prob=[0.1;0.90];
    otherwise
end
        


% if the postsynaptic cell has good connection, switch to 
disconnected_threshold_good = 0.2;connected_threshold_good = 0.5;
quantiles_prob_good=[0.1;0.90];

multiplier_connected=2;

% specify the parameters:
% som_params
% regular_params
% disconnected_params

% determine how many batches are used in fitting the model, default is zero 
num_trace_back = 2;

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



%% Online design:
model_type=1; % 1 working, 2 ss 3 fuparamell
design_type_multi=1; %0 random 1 optimal 
design_type_single=1; %0 random 1 optimal 
synchrony_test=1;
small_psc_events=0;
post_cell_type=1;

disconnected_indicators=ones(n_cell_this_plane,1);
% if design_type == 1 % optimal design
%     disp('Use optimal design')
% elseif design_type==0 % random design
%     disp('Use all random design')
% end
switch model_type
    case 1% working model
        disp('working model')
        spike_indicator=false;
    case 2%
        disp('working model with spike & slab')
        spike_indicator=true;
    case 3 %
        disp(' full model with first event time')
    otherwise     
end

switch assignment_type
    case 1
        disp('gamma only')
    case 2 %
        disp('gamma and spike probability')
    case 3 %
        disp('gamma and gain variance')
end
% tuning parameters: model_type, design_type, assignment_type
run('./Simulation_experiment.m')
%%
%gamma_related=gamma_truth(related_cell_list);
%gain_related=gain_truth(related_cell_list);


    save(strcat('./matfiles/Sep25/','Sim', num2str(i_sim),'Seed',num2str(i_seed),'.mat'),...
        'variational_params_path','gamma_path','gain_path','var_gamma_path',...
        'mpp_connected', 'trials_locations_connected','trials_powers_connected',...
        'mpp_disconnected', 'trials_locations_disconnected','trials_powers_disconnected',...
        'mpp_undefined', 'trials_locations_undefined','trials_powers_undefined',...
        'undefined_cells', 'potentially_disconnected_cells', 'potentially_connected_cells',...
    'dead_cells', 'alive_cells')
