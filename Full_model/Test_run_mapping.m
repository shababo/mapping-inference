addpath(genpath('../../mapping-inference'),genpath('../../odessa-beta-beta'));
%% Run get_experiment_setup:

disp('Test get_experiment_setup')
experiment_setup=get_experiment_setup();
experiment_setup.is_exp = 0;
experiment_setup.enable_user_breaks = 0;
%% Generate cells
disp('Test generate_neurons')
experiment_setup.sim=get_simulation_setup();
experiment_setup.neurons=generate_neurons(experiment_setup);
%% Initialize neighbourhoods 
disp('Test create_neighbourhoods')
neighbourhoods = create_neighbourhoods(experiment_setup);
    
%% 
disp('Test data collection')
