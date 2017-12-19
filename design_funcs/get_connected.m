
function connected_profile= get_connected()

% basics
connected_profile=struct;
connected_profile.group_ID='connected';

connected_profile.design_function=@design_connected;

connected_profile.design_func_params=struct;
connected_profile.design_func_params.candidate_grid_params=struct;
connected_profile.design_func_params.candidate_grid_params.radius=[5 12.5];
connected_profile.design_func_params.candidate_grid_params.number=[12 16];
connected_profile.design_func_params.candidate_grid_params.grid_type='ring'; % 2d ring or 3d sphere


connected_profile.design_func_params.trials_params=struct;
connected_profile.design_func_params.trials_params.replicates=1;
connected_profile.design_func_params.trials_params.spots_per_trial=1;
connected_profile.design_func_params.trials_per_cell=20;
connected_profile.design_func_params.trials_params.power_levels=30:10:100;
connected_profile.design_func_params.trials_params.stim_design='Optimal';
connected_profile.design_func_params.trials_params.MCsamples_for_posterior=50;
connected_profile.design_func_params.trials_params.trials_per_batch=400;   
connected_profile.design_func_params.trials_params.num_stim_sites=4;   
connected_profile.design_func_params.trials_params.min_gap_stim=9.5;   
connected_profile.design_func_params.trials_params.bounds.gain=[0.005 0.03];

% Random, Nuclei, or Optimal
%   whether to conduct more trials on low PR cells
connected_profile.design_func_params.trials_params.weighted_indicator=true;


connected_profile.psc_detect_function = @run_oasis;
connected_profile.inference_function = @inference_connected;

connected_profile.inference_params=struct;
connected_profile.inference_params.likelihood=@lif_glm_firstspike_loglikelihood_for_VI;
connected_profile.inference_params.maxit=1e4;
connected_profile.inference_params.MCsamples_for_gradient=50;
connected_profile.inference_params.convergence_threshold=1e-3;
connected_profile.inference_params.step_size=1;
connected_profile.inference_params.step_size_max=2;
connected_profile.inference_params.MCsamples_for_posterior=50;
connected_profile.inference_params.recent_batches=2;
connected_profile.inference_params.bounds=struct;
connected_profile.inference_params.bounds.PR=[0.05 1];
connected_profile.inference_params.bounds.gain=[0.005 0.03];
connected_profile.inference_params.bounds.spike_indicator=false;


connected_profile.regroup_function=struct;
connected_profile.regroup_function.disconnected=@connected_to_disconnected;
connected_profile.regroup_function.alive=@connected_to_alive;
connected_profile.regroup_func_params=struct;
connected_profile.regroup_func_params.connected_threshold=0.5;
connected_profile.regroup_func_params.disconnected_threshold=0.2;
connected_profile.regroup_func_params.quantile_prob=[0.05 0.95];
connected_profile.regroup_func_params.regroup_type='Quantiles'; % Quantiles or NonzeroProb
connected_profile.regroup_func_params.singlespot_threshold=0.2;% certain proportion of cells 
connected_profile.regroup_func_params.change_threshold =0.05;




