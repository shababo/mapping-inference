function secondary_profile= get_secondary()

% basics
secondary_profile=struct;
secondary_profile.group_ID='secondary';

secondary_profile.inference_params=struct;
secondary_profile.inference_params.likelihood=@calculate_loglikelihood_bernoulli;
secondary_profile.inference_params.maxit=1e4;
secondary_profile.inference_params.MCsamples_for_gradient=50;
secondary_profile.inference_params.convergence_threshold=1e-2;
secondary_profile.inference_params.step_size=1;
secondary_profile.inference_params.step_size_max=2;
secondary_profile.inference_params.MCsamples_for_posterior=50;
secondary_profile.inference_params.recent_batches=2;
secondary_profile.inference_params.bounds=struct;
secondary_profile.inference_params.bounds.PR=[0.05 1];
secondary_profile.inference_params.bounds.gain=[0.005 0.03];
secondary_profile.inference_params.bounds.spike_indicator=false;


secondary_profile.regroup_func_params=struct;
secondary_profile.regroup_func_params.connected_threshold=0.5;
secondary_profile.regroup_func_params.disconnected_threshold=0.2;
secondary_profile.regroup_func_params.quantile_prob=[0.1 0.9];
secondary_profile.regroup_func_params.regroup_type='Quantiles'; % Quantiles or NonzeroProb
secondary_profile.regroup_func_params.singlespot_threshold=0.2;% certain proportion of cells 
secondary_profile.regroup_func_params.change_threshold =0.05;




