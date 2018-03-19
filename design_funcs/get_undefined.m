function undefined_profile= get_undefined()

% basics
undefined_profile=struct;
undefined_profile.group_ID='undefined';

undefined_profile.design_function=@design_undefined;

undefined_profile.design_func_params=struct;
undefined_profile.design_func_params.candidate_grid_params=struct;
undefined_profile.design_func_params.candidate_grid_params.radius=[2.5 6.5];
undefined_profile.design_func_params.candidate_grid_params.number=[8 12];
undefined_profile.design_func_params.candidate_grid_params.grid_type='ring'; % 2d ring or 3d sphere


undefined_profile.design_func_params.trials_params=struct;
undefined_profile.design_func_params.trials_params.replicates=1;
undefined_profile.design_func_params.trials_params.spots_per_trial=3;
undefined_profile.design_func_params.trials_per_cell=10;
undefined_profile.design_func_params.trials_params.power_levels=50:10:100;
undefined_profile.design_func_params.trials_params.stim_design='Optimal';
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
% undefined_profile.inference_params.likelihood=@lif_glm_firstspike_loglikelihood_for_VI;
undefined_profile.inference_params.likelihood=@lif_glm_firstevent_loglikelihood_for_VI;
undefined_profile.inference_params.maxit=500;
undefined_profile.inference_params.MCsamples_for_gradient=100;
undefined_profile.inference_params.convergence_threshold=1e-3;
undefined_profile.inference_params.step_size=1;
undefined_profile.inference_params.step_size_max=1;
undefined_profile.inference_params.MCsamples_for_posterior=50;
undefined_profile.inference_params.recent_batches=2;
undefined_profile.inference_params.bounds=struct;
undefined_profile.inference_params.bounds.PR=[0.05 1];
undefined_profile.inference_params.bounds.gain=[0.001 0.1];
undefined_profile.inference_params.bounds.spike_indicator=false;
undefined_profile.inference_params.bounds.delay_mu=[0 60];
undefined_profile.inference_params.bounds.delay_sigma=[0.1 10];


undefined_profile.regroup_function=struct;
undefined_profile.regroup_function.connected=@undefined_to_connected;
undefined_profile.regroup_function.disconnected=@undefined_to_disconnected;

undefined_profile.regroup_func_params=struct;
undefined_profile.regroup_func_params.connected_threshold=0.5;
undefined_profile.regroup_func_params.disconnected_threshold=0.2;
undefined_profile.regroup_func_params.quantile_prob=[0.05 0.95];
undefined_profile.regroup_func_params.regroup_type='Quantiles'; % Quantiles or NonzeroProb
undefined_profile.regroup_func_params.singlespot_threshold=25;%0.2;% certain proportion of cells 
undefined_profile.regroup_func_params.change_threshold =0.05;




