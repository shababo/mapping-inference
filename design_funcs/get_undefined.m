function undefined_profile= get_undefined()

% basics
undefined_profile=struct;
undefined_profile.group_ID='undefined';

undefined_profile.design_function=@design_undefined;

undefined_profile.design_func_params=struct;
undefined_profile.design_func_params.candidate_grid_params=struct;
undefined_profile.design_func_params.candidate_grid_params.max_radius=10; % for each dim
undefined_profile.design_func_params.candidate_grid_params.number=[1 5 10];

% undefined_profile.design_func_params.candidate_grid_params.number=[];
undefined_profile.design_func_params.candidate_grid_params.grid_type='ring'; % ring, or random


undefined_profile.design_func_params.trials_params=struct;
undefined_profile.design_func_params.trials_params.num_stim_sites=4; % for each dim
undefined_profile.design_func_params.trials_params.replicates=1;
undefined_profile.design_func_params.trials_params.spots_per_trial=1;
undefined_profile.design_func_params.trials_per_cell=10;
undefined_profile.design_func_params.trials_params.power_levels=20:5:50;
undefined_profile.design_func_params.trials_params.stim_design='Random';
undefined_profile.design_func_params.trials_params.MCsamples_for_posterior=50;
undefined_profile.design_func_params.trials_params.trials_per_batch=300;   


% Random, Nuclei, or Optimal
%   whether to conduct more trials on low PR cells
undefined_profile.design_func_params.trials_params.weighted_indicator=true;

undefined_profile.psc_detect_function = @run_oasis;
undefined_profile.inference_function = @inference_undefined;

undefined_profile.inference_params=struct;
% Three likelihoods: 
% undefined_profile.inference_params.likelihood=@calculate_loglikelihood_bernoulli;
undefined_profile.inference_params.likelihood=@lif_glm_firstspike_loglikelihood_for_VI;
undefined_profile.inference_params.maxit=500;
undefined_profile.inference_params.MCsamples_for_gradient=40;
undefined_profile.inference_params.convergence_threshold=2e-3;
undefined_profile.inference_params.step_size=1;
undefined_profile.inference_params.step_size_max=0.1;
undefined_profile.inference_params.MCsamples_for_posterior=50;
undefined_profile.inference_params.recent_batches=2;
undefined_profile.inference_params.event_range=[1 300];
undefined_profile.inference_params.bounds=struct;

undefined_profile.regroup_function=struct;
undefined_profile.regroup_function.connected=@undefined_to_connected;
undefined_profile.regroup_function.disconnected=@undefined_to_disconnected;

undefined_profile.regroup_func_params=struct;
undefined_profile.regroup_func_params.connected_threshold=0.5;
undefined_profile.regroup_func_params.disconnected_threshold=0.2;
undefined_profile.regroup_func_params.quantile_prob=[0.05 0.5 0.95];
undefined_profile.regroup_func_params.regroup_type='Quantiles'; % Quantiles or NonzeroProb
undefined_profile.regroup_func_params.singlespot_threshold=0.2;%0.2;% certain proportion of cells 
undefined_profile.regroup_func_params.change_threshold =0.05;
undefined_profile.regroup_func_params.undefined_threshold=6;



