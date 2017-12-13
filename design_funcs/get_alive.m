function alive_profile= get_alive()

% basics
alive_profile=struct;
alive_profile.group_ID='alive';

alive_profile.design_function=@design_alive;

connected_profile.design_func_params=struct;
connected_profile.design_func_params.candidate_grid_params=struct;
connected_profile.design_func_params.candidate_grid_params.radius=[];
connected_profile.design_func_params.candidate_grid_params.number=[0];
connected_profile.design_func_params.candidate_grid_params.grid_type='ring'; % 2d ring or 3d sphere

connected_profile.design_func_params.trials_params=struct;
connected_profile.design_func_params.trials_params.replicates=1;
connected_profile.design_func_params.trials_params.spots_per_trial=1;
connected_profile.design_func_params.trials_per_cell=20;
connected_profile.design_func_params.trials_params.power_levels=30:10:100;
connected_profile.design_func_params.trials_params.stim_design='Optimal';
connected_profile.design_func_params.trials_params.MCsamples_for_posterior=50;
connected_profile.design_func_params.trials_params.trials_per_batch=200;   

%   whether to conduct more trials on low PR cells
connected_profile.design_func_params.trials_params.weighted_indicator=true;

connected_profile.psc_detect_function = @run_oasis;
connected_profile.inference_function = @no_inference;