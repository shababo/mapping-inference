% function [this_neighbourhood]=inference_undefined_dev(...
%     experiment_query_this_group,this_neighbourhood,group_profile,experiment_setup)
% 
group_ID='undefined';
% [cells_this_group, group_cell_IDs] = get_group_inds(this_neighbourhood,group_ID);
% cells_this_group = find(cells_this_group);
% number_cells_this_group=length(cells_this_group);
% number_cells_all= length(this_neighbourhood.neurons);
prior_info=experiment_setup.prior_info;
group_profile=experiment_setup.groups.(group_ID);
% 
% %indicators_remained = find(ismember([mpp_undefined(:).batch],iter-(0:num_trace_back) ));
% number_of_trials = length(experiment_query_this_group.trials);

% mpp_remained=extract_mpp(experiment_query_this_group.trials);
% designs_remained = get_stim_size(group_ID,experiment_query_this_group.trials,this_neighbourhood);

mpp_remained=extract_mpp(available_trials(cell_with_trials(i)).trials);
designs_remained = [available_trials(cell_with_trials(i)).trials(:).stimulation]';

%
% include all cells that have been stimulated:
stim_threshold=prior_info.induced_intensity.minimum_stim_threshold/group_profile.inference_params.bounds.gain(2);
i_stimulated_cell_list= find(sum(designs_remained>stim_threshold,1)>0);
number_of_stim_cells=length(i_stimulated_cell_list);
designs_remained=designs_remained(:,i_stimulated_cell_list);

%%
% Update variational and prior distribution


% batch_ID=this_neighbourhood.batch_ID;
% neurons=this_neighbourhood.neurons(i_stimulated_cell_list);
% properties={'PR_params','gain_params'};summary_stat={'pi_logit','alpha','beta'};
% temp_output=grab_values_from_neurons(Inf,neurons,properties,summary_stat);

variational_params=struct;
% There are two parameters in the mean model: a and phi
% There are two new parameter in the variance model: c and d;
% And gamma

% We need a new data structure for saving this parameter!

% Log-normal for a,b,c,d
% Logit-normal for gain and PR (due to upper bound)
% Specify their parameters as well as their distibutions:
variational_params(1).a = struct;
variational_params(1).a(1).params = [-1 1];
variational_params(1).a(1).family = 'log-normal';

variational_params(1).b = struct;
variational_params(1).b(1).params = [-1 1];
variational_params(1).b(1).family = 'log-normal';

variational_params(1).c = struct;
variational_params(1).c(1).params = [-1 1];
variational_params(1).c(1).family = 'log-normal';

variational_params(1).d = struct;
variational_params(1).d(1).params = [-1 1];
variational_params(1).d(1).family = 'log-normal';

variational_params(1).gain = struct;
variational_params(1).gain(1).params = [0 1];
variational_params(1).gain(1).family = 'logit-normal';

variational_params(1).PR = struct;
variational_params(1).PR(1).params = [0 1];
variational_params(1).PR(1).family = 'logit-normal';

prior_params=variational_params;


%%
%prior_params=variational_params_path(max(iter-num_trace_back,1),cell_list);


[parameter_history] = fit_VI_dev(...
    designs_remained, mpp_remained, experiment_setup.patched_neuron.background_rate,...
    variational_params,prior_params,...
    group_profile.inference_params,prior_info);


    % need a function that writes the values into the parameter history in
    % this_neighbourhood 



for i_cell = 1:number_of_stim_cells
    
    current_params=reformat_to_neurons(parameter_history(end,i_cell),'gamma','spiked_logit_normal');
    
%     group_profile=experiment_setup.groups.(this_neighbourhood.neurons(i_cell).group_ID);
    bounds= group_profile.inference_params.bounds.PR;
    quantile_prob=group_profile.regroup_func_params.quantile_prob;
    this_neighbourhood.neurons(i_stimulated_cell_list(i_cell)).PR_params(batch_ID)=calculate_posterior(...
        current_params,bounds,quantile_prob);
   
    current_params=reformat_to_neurons(parameter_history(end,i_cell),'gain','spiked_logit_normal');
    bounds= group_profile.inference_params.bounds.gain;
    quantile_prob=group_profile.regroup_func_params.quantile_prob;
    this_neighbourhood.neurons(i_stimulated_cell_list(i_cell)).gain_params(batch_ID)=calculate_posterior(...
        current_params,bounds,quantile_prob);
    
end

                
end
