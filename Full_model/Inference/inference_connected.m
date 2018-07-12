function [neighbourhood]=inference_connected(...
    experiment_query_this_group,neighbourhood,group_profile,experiment_setup)

%  experiment_query_this_group=experiment_query.(this_group);
%  this_neighbourhood= neighbourhood;
%  experiment_query.(this_group),neighbourhood,group_profile, experiment_setup
%%

group_ID=group_profile.group_ID;
[cells_this_group, group_cell_IDs] = get_group_inds(neighbourhood,group_ID);
cells_this_group = find(cells_this_group);
number_cells_this_group=length(cells_this_group);
number_cells_all= length(neighbourhood.neurons);
prior_info=experiment_setup.prior_info;

inference_params=group_profile.inference_params;
%indicators_remained = find(ismember([mpp_undefined(:).batch],iter-(0:num_trace_back) ));
number_of_trials = length(experiment_query_this_group.trials);

% Need a function that graph the mpp from the experiment_query
% note: this_experiment_query contains the group information 
trials = experiment_query_this_group.trials;
neurons=neighbourhood.neurons;
[variational_params, prior_params]=initialize_params_VI(neurons);
[variational_params]=initialize_PR_VI(variational_params,neurons,trials,prior_info);
%%

% prior_info.prior_parameters.boundary_params= [30 30 70];

[parameter_history,elbo_rec,loglklh_rec] = fit_VI(...
      trials,neurons, experiment_setup.patched_neuron.background_rate, ...
    variational_params,prior_params,...
    inference_params,prior_info);

%% Update the parameters in neighbourhood 
   
batch_ID=neighbourhood.batch_ID;
quantile_prob=group_profile.regroup_func_params.quantile_prob;
       
for i_cell = 1:number_cells_all
    
    neighbourhood.neurons(i_cell).params(batch_ID)=parameter_history(end,i_cell);
    neighbourhood.neurons(i_cell).posterior_stat(batch_ID)=...
        calculate_posterior(parameter_history(end,i_cell),quantile_prob);
        
end
%%

% [clusters_of_cells] = find_clusters(stim_all, 1:num_cells_nhood, stim_threshold);

