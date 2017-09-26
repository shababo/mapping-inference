addpath(genpath('../../../mapping-inference'));
%% Specify the setting in this simulation 

for i_sim = 0:4
    grid_type_global = 1; % grid for connected cells, 1: circles; 2: lines;
    multi_stim_type_global = 1; % sim spots for multi-spot stimulation, 1: targeted; 2: random
    single_spot_threshold_global=9; % number of cells before switching to single cell stimulation 
    power_type_global = 1; % power rules, 1: fixed; 2: random 
    if i_sim == 0 % default 
    elseif i_sim == 1 % line v.s. circle 
        grid_type_global=2;
    elseif i_sim == 2
        power_type_global =2;
    elseif i_sim == 3
        single_spot_threshold_global = 15;
    elseif i_sim == 4
        multi_stim_type_global = 2;
    end
        
    for i_seed = 1:20
        rng(i_seed,'twister');
       


%% Generate cellular parameters
run('./Simulation_parameters.m')
%% Preprocessing

grid_type = grid_type_global; % 1: circles; 2: lines;
multi_stim_type=multi_stim_type_global; % 1: target; 2: random;
run('./Simulation_preprocessing.m')

% Sim 1: three optimal design for connected cells 
% Sim 2: nothing
% Sim 3: nothing
% Sim 4: targeted cell-killing v.s. random stim locations 
%% Designing experiment
%-------------------------------------------------%

% Online design:
%% Parameters in the design stage
% Sim 1: nothing (changes are made in location selection)
% Sim 2: Enable random power level for connected cells (in random design)
% Sim 3: New assignment rules, and explore switch to single-spot stim
% Sim 4: nothing 

run('./Simulation_design_choices.m')

n_replicates=1; % conduct two replicates for each trial
K_undefined=5; % each cell appears approximately 10*2 times
K_disconnected=5; % each cell appears approximately 10*2 times
K_connected=10; % each cell appears approximately 10*2 times

single_spot_threshold=single_spot_threshold_global; % switch to single spot stimulation (this can be a function of n_spots_per_trial
trial_max=2000;

% threshold for the group movement:
disconnected_threshold = 0.2;disconnected_confirm_threshold = 0.2;
connected_threshold = 0.5;connected_confirm_threshold = 0.5;
change_threshold=0.05; % for potentially connected cells (since the estimated gamma are important in this case)

new_rule=false;
%% Online design:
% Sim 1: nothing (changes are made in location selection)
% Sim 2: Enable random power level for connected cells (in random design)
% Sim 3: New assignment rules 
% Sim 4: nothing 
power_type = power_type_global;
run('./Simulation_exp.m')
%% Outputs:
    
    save(strcat('./matfiles/Sep25/','Sim', num2str(i_sim),'Seed',num2str(i_seed),'.mat'),...
        'variational_params_path','gamma_path','var_gamma_path',...
        'mpp_connected', 'trials_locations_connected','trials_powers_connected',...
        'mpp_disconnected', 'trials_locations_disconnected','trials_powers_disconnected',...
        'mpp_undefined', 'trials_locations_undefined','trials_powers_undefined',...
        'undefined_cells', 'potentially_disconnected_cells', 'potentially_connected_cells',...
    'dead_cells', 'alive_cells')
        
        
    end
end