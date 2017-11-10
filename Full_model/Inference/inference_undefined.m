function [this_neighbourhood]=inference_undefined(...
    this_experiment_query,this_neighbourhood,group_profile,prior_info)

group_type_ID=group_profile.group_type_ID;
cells_this_group=index([this_neighbourhood.neurons(:).group_type_ID]==group_type_ID);
number_cells_this_group=length(cells_this_group);
number_cells_all= length(this_neighbourhood.neurons);


%indicators_remained = find(ismember([mpp_undefined(:).batch],iter-(0:num_trace_back) ));
number_of_trials = length(this_experiment_query.trials);

% Need a function that graph the mpp from the experiment_query
% note: this_experiment_query contains the group information 
mpp_remained=extract_mpp(this_experiment_query);
designs_remained = get_stim_size(this_experiment_query.trials,this_neighbourhood);


% include all cells that have been stimulated:
stimulated_cell_list= find(sum(designs_remained>prior_info.induced_intensity.stim_threshold,1)>0);
designs_remained=designs_remained(:,stimulated_cell_list);
% Update variational and prior distribution

variational_params=grab_recent_value(this_neighbourhood.neurons(stimulated_cell_list));
prior_params=grab_recent_value(this_neighbourhood.neurons(stimulated_cell_list));

%prior_params=variational_params_path(max(iter-num_trace_back,1),cell_list);


lklh_func=group_profile.inference_params.likelihood

[parameter_history] = fit_VI(...
    designs_remained, mpp_remained, variational_params,prior_params,
    group_profile.inference_params);

    % need a function that writes the values into the parameter history in
    % this_neighbourhood 
[parameter_temp] = calculate_posterior(parameter_history(end,:),gamma_bound,gain_bound,quantiles_prob,spike_indicator);
variational_params_path(iter+1,cell_list)=parameter_history(end,:);

end
