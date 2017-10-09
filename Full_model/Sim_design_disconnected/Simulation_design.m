addpath(genpath('../../../mapping-inference'));
%% Specify the setting in this simulation 

for i_sim = [12 15]
    grid_type_global = 1; % grid for connected cells, 1: circles; 2: lines;
    multi_stim_type_global = 1; % sim spots for multi-spot stimulation, 1: targeted; 2: random
    single_spot_threshold_global=9; % number of cells before switching to single cell stimulation 
    power_type_global = 1; % power rules, 1: fixed; 2: random 
    gain_type_global = 1; % 1 default; 2: long-tail in low gains; 
    multi_power_type_global=1; % 1: default 50mV; 2: high power 150mV;
    
good_prior_global=1; % 1: good prior; 2: uninformative prior
random_power_for_gain_global=1; % 1: use random power; 2: use chosen power;
fit_gain_global=1; % 1: fit gains in working model; 2: don't fit gain.  

    if i_sim == 0 % default 
    elseif i_sim == 1 % line v.s. circle 
        grid_type_global=2;
    elseif i_sim == 2
        power_type_global =2;
    elseif i_sim == 3
        single_spot_threshold_global = 15;
    elseif i_sim == 4
        multi_stim_type_global = 2;
    elseif i_sim == 5
        gain_type_global=2;power_type_global =2; % long-tail weakly expressed & random power levels for connected cells
    elseif i_sim == 6
        gain_type_global=2;power_type_global =2;
        multi_power_type_global=2; 
    elseif i_sim==10 % bad data, good prior, old method
        gain_type_global=2; fit_gain_global=2;
        power_type_global=2;good_prior_global=1; 
        
    elseif i_sim==11 % bad data, good prior, new method (random)
        gain_type_global=2; fit_gain_global=1;
        power_type_global=2;good_prior_global=1; 
        
    elseif i_sim==12 % bad data, good prior, new method (selected)
        gain_type_global=2; fit_gain_global=1;
        power_type_global=2;good_prior_global=1; 
        random_power_for_gain_global=2;
        
    elseif i_sim==13 % bad data, bad prior, old method
        gain_type_global=2; fit_gain_global=2;
        power_type_global=2;good_prior_global=2;
        
    elseif i_sim==14 % bad data, bad prior, new method (random)
        gain_type_global=2; fit_gain_global=1;
        power_type_global=2;good_prior_global=2;
        
    elseif i_sim==15 % bad data, bad prior, new method (selected)
        gain_type_global=2; fit_gain_global=1;
        power_type_global=2;good_prior_global=2;
        random_power_for_gain_global=2;
    end
       
    for i_seed = 1:20
        rng(i_seed,'twister');
       
%% Generate cellular parameters
gain_type =gain_type_global;
good_prior=good_prior_global;
run('./Simulation_parameters.m')
%% Preprocessing
%gain_truth(:)=0.005;
fit_gain=fit_gain_global;
grid_type = grid_type_global; % 1: circles; 2: lines;
multi_stim_type=multi_stim_type_global; % 1: target; 2: random;
multi_power_type=multi_power_type_global;
if fit_gain==1
run('./Simulation_preprocessing_gain.m')
else
run('./Simulation_preprocessing.m')
    
end
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
if fit_gain==1
run('./Simulation_design_choices_gain.m')
else
run('./Simulation_design_choices.m')
end    
n_replicates=1; % conduct two replicates for each trial
K_undefined=20; % each cell appears approximately 10*2 times
K_disconnected=5; % each cell appears approximately 10*2 times
K_connected=5; % each cell appears approximately 10*2 times
n_spots_per_trial = 4;
single_spot_threshold=6; % switch to single spot stimulation (this can be a function of n_spots_per_trial
trial_max=2000;

% threshold for the group movement:
disconnected_threshold = 0.2;disconnected_confirm_threshold = 0.2;
connected_threshold = 0.5;connected_confirm_threshold = 0.5;
change_threshold=0.1; % for potentially connected cells (since the estimated gamma are important in this case)

new_rule=false;
%% Online design:
% Sim 1: nothing (changes are made in location selection)
% Sim 2: Enable random power level for connected cells (in random design)
% Sim 3: New assignment rules 
% Sim 4: nothing 
power_type = power_type_global;
random_power_for_gain=random_power_for_gain_global;

if fit_gain==1
    run('./Simulation_exp_gain.m')
else
    run('./Simulation_exp.m')
end
%
% find(alive_cells{iter})
% find(gamma_related>0)
%
% gain_related

%
% find(gamma_truth(target_cell_list.primary)>0)
% temp_gain=gain_truth(target_cell_list.primary);
% temp_gain(find(gamma_truth(target_cell_list.primary)>0))
% find(alive_cells{iter})
% Outputs:
    %%
    save(strcat('./matfiles/Sep25/','Sim', num2str(i_sim),'Seed',num2str(i_seed),'.mat'),...
        'variational_params_path','gamma_path','gain_path','var_gamma_path',...
        'mpp_connected', 'trials_locations_connected','trials_powers_connected',...
        'mpp_disconnected', 'trials_locations_disconnected','trials_powers_disconnected',...
        'mpp_undefined', 'trials_locations_undefined','trials_powers_undefined',...
        'undefined_cells', 'potentially_disconnected_cells', 'potentially_connected_cells',...
    'dead_cells', 'alive_cells')
        
        
    end
end