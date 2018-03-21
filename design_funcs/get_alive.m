function alive_profile= get_alive()

% basics
alive_profile=struct;
alive_profile.group_ID='alive';

alive_profile.design_function=@design_alive;

alive_profile.design_func_params=struct;
alive_profile.design_func_params.candidate_grid_params=struct;
alive_profile.design_func_params.candidate_grid_params.radius=[];
alive_profile.design_func_params.candidate_grid_params.number=[0];
alive_profile.design_func_params.candidate_grid_params.grid_type='ring'; % 2d ring or 3d sphere

alive_profile.design_func_params.trials_params=struct;
alive_profile.design_func_params.trials_params.replicates=1;
alive_profile.design_func_params.trials_params.spots_per_trial=1;
alive_profile.design_func_params.trials_per_cell=20;
alive_profile.design_func_params.trials_params.power_levels=30:10:100;
alive_profile.design_func_params.trials_params.stim_design='Optimal';
alive_profile.design_func_params.trials_params.MCsamples_for_posterior=50;
alive_profile.design_func_params.trials_params.trials_per_batch=200;   
alive_profile.design_func_params.trials_params.bounds.gain=[0.005 0.03];
alive_profile.design_func_params.trials_params.num_stim_sites = 1;
alive_profile.design_func_params.trials_params.min_gap_stim = 0;


%   whether to conduct more trials on low PR cells
alive_profile.design_func_params.trials_params.weighted_indicator=true;

alive_profile.psc_detect_function = @run_oasis;
alive_profile.inference_function = @no_inference;