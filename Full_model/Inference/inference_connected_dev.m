function [this_neighbourhood]=inference_connected_dev(...
    experiment_query_this_group,this_neighbourhood,group_profile,experiment_setup)

%  experiment_query_this_group=experiment_query.(this_group);
%  this_neighbourhood= neighbourhood;
%  experiment_query.(this_group),neighbourhood,group_profile, experiment_setup

group_ID=group_profile.group_ID;
i_cell_group_to_nhood= find(get_group_inds(this_neighbourhood,group_ID));
number_cells_this_group=length(i_cell_group_to_nhood);
num_cells_nhood= length(this_neighbourhood.neurons);
prior_info=experiment_setup.prior_info;

%indicators_remained = find(ismember([mpp_undefined(:).batch],iter-(0:num_trace_back) ));
number_of_trials = length(experiment_query_this_group.trials);

% Need a function that graph the mpp from the experiment_query
% note: this_experiment_query contains the group information 
mpp_all=extract_mpp(experiment_query_this_group.trials);
stim_all = get_stim_size(group_ID,experiment_query_this_group.trials,this_neighbourhood);


% include all cells that have been stimulated:
stim_threshold=prior_info.induced_intensity.minimum_stim_threshold/group_profile.inference_params.bounds.gain(2);
%stimulated_cell_list= find(sum(stim_all>stim_threshold,1)>0);

%number_of_stim_cells=length(stimulated_cell_list);
% designs_remained=designs_remained(:,stimulated_cell_list);

[clusters_of_cells] = find_clusters(stim_all, 1:num_cells_nhood, stim_threshold);


% Update variational and prior distribution
batch_ID=this_neighbourhood.batch_ID;
neurons=this_neighbourhood.neurons(:);
properties={'PR_params'};summary_stat={'pi_logit','alpha','beta'};
temp_output=grab_values_from_neurons(Inf,neurons,properties,summary_stat);

variational_params=struct;
temp=num2cell(temp_output.PR_params.pi_logit); [variational_params(1:num_cells_nhood).p_logit]=temp{:};
temp=num2cell(temp_output.PR_params.alpha); [variational_params(1:num_cells_nhood).alpha]=temp{:};
temp=num2cell(temp_output.PR_params.beta); [variational_params(1:num_cells_nhood).beta]=temp{:};

properties={'gain_params','delay_mu_params','delay_sigma_params'};
summary_stat={'alpha','beta'};
temp_output=grab_values_from_neurons(Inf,neurons,properties,summary_stat);
temp=num2cell(temp_output.gain_params.alpha); [variational_params(1:num_cells_nhood).alpha_gain]=temp{:};
temp=num2cell(temp_output.gain_params.beta); [variational_params(1:num_cells_nhood).beta_gain]=temp{:};
temp=num2cell(temp_output.delay_mu_params.alpha); [variational_params(1:num_cells_nhood).alpha_m]=temp{:};
temp=num2cell(temp_output.delay_mu_params.beta); [variational_params(1:num_cells_nhood).beta_m]=temp{:};
temp=num2cell(temp_output.delay_sigma_params.alpha); [variational_params(1:num_cells_nhood).alpha_s]=temp{:};
temp=num2cell(temp_output.delay_sigma_params.beta); [variational_params(1:num_cells_nhood).beta_s]=temp{:};
% prior_params=variational_params;

%prior_params=variational_params_path(max(iter-num_trace_back,1),cell_list);

for i_cluster= 1:length(clusters_of_cells)
    %
    % find the trials that are relevant to thiscell
    active_trials=find(sum(stim_all(:,(clusters_of_cells{i_cluster})),2)>stim_threshold);
    neighbour_list= find(sum(stim_all(active_trials,:)>stim_threshold,1)>0);
    if ~isempty(neighbour_list)
        
        var_params=variational_params(neighbour_list);
        prior_params=var_params;
        stim_size=stim_all(active_trials,neighbour_list);
        mpp=mpp_all(active_trials);
        
        
        [parameter_history,loglklh_rec] = fit_VI_dev3(...
            stim_size, mpp, experiment_setup.patched_neuron.background_rate,...
            var_params,prior_params,...
            group_profile.inference_params,prior_info);

       
         
        for i_cell = 1:length(neighbour_list)
            current_params=reformat_to_neurons(parameter_history(end,i_cell),'gamma','spiked_logit_normal');
%             group_profile=experiment_setup.groups.(this_neighbourhood.neurons(i_cell).group_ID);
            bounds= group_profile.inference_params.bounds.PR;
            quantile_prob=group_profile.regroup_func_params.quantile_prob;
            this_neighbourhood.neurons(neighbour_list(i_cell)).PR_params(batch_ID)=calculate_posterior(...
                current_params,bounds,quantile_prob);
            current_params=reformat_to_neurons(parameter_history(end,i_cell),'gain','spiked_logit_normal');
%             group_profile=experiment_setup.groups.(this_neighbourhood.neurons(i_cell).group_ID);
            bounds= group_profile.inference_params.bounds.gain;
            quantile_prob=group_profile.regroup_func_params.quantile_prob;
            this_neighbourhood.neurons(neighbour_list(i_cell)).gain_params(batch_ID)=calculate_posterior(...
                current_params,bounds,quantile_prob);
            
            current_params=reformat_to_neurons(parameter_history(end, this_cell),'delay_mu','spiked_logit_normal');
            bounds= [delay_mu_bound.low delay_mu_bound.up];
            this_neighbourhood.neurons(i_stimulated_cell_list(i_cell)).delay_mu_params=calculate_posterior(...
                current_params,bounds,quantile_prob);
            
            current_params=reformat_to_neurons(parameter_history(end, this_cell),'delay_sigma','spiked_logit_normal');
            bounds= [delay_sigma_bound.low delay_sigma_bound.up];
            this_neighbourhood.neurons(i_stimulated_cell_list(i_cell)).delay_sigma_params=calculate_posterior(...
                current_params,bounds,quantile_prob);
            
        end

        
        
    end
end




                
end
