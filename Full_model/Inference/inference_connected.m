function [neighbourhood]=inference_connected(...
    experiment_query_this_group,neighbourhood,group_profile,experiment_setup)
%%
%  experiment_query_this_group=experiment_query.(this_group);
%  neighbourhood= neighbourhoods(1);
%  experiment_query.(this_group),neighbourhood,group_profile, experiment_setup
%%
group_ID=group_profile.group_ID;
[cells_this_group, group_cell_IDs] = get_group_inds(neighbourhood,group_ID);
cells_this_group = find(cells_this_group);
number_cells_this_group=length(cells_this_group);
number_cells_all= length(neighbourhood.neurons);
prior_info=experiment_setup.prior_info;
background_rate= experiment_setup.patched_neuron.background_rate;

inference_params=group_profile.inference_params;
%indicators_remained = find(ismember([mpp_undefined(:).batch],iter-(0:num_trace_back) ));
number_of_trials = length(experiment_query_this_group.trials);

% Need a function that graph the mpp from the experiment_query
% note: this_experiment_query contains the group information 
trials = experiment_query_this_group.trials;
neurons=neighbourhood.neurons;
[variational_params, prior_params,trials]=initialize_params_VI(neurons,trials,prior_info); 
% [variational_params]=initialize_PR_VI(variational_params,neurons,trials,prior_info,inference_params,background_rate);
% variational_params.shapes.dist='logit-normal';
% variational_params.delay_mu.mean=-6;variational_params.delay_mu.log_sigma=-4;
% variational_params.delay_sigma.mean=-6;variational_params.delay_sigma.log_sigma=-4;
% prior_params=variational_params;
% variational_params.gain.mean = log( (0.0159-0.005)/(0.06-0.0159));
%%
% prior_info.prior_parameters.boundary_params= [30 30 70];
% inference_params.maxit=500;
% inference_params.step_size_max=0.1;

[parameter_history,elbo_rec] = fit_VI(...
      trials,neurons, experiment_setup.patched_neuron.background_rate, ...
    variational_params,prior_params,...
    inference_params,prior_info);
%%

batch_ID=neighbourhood.batch_ID;
quantile_prob=group_profile.regroup_func_params.quantile_prob;
for i_cell = 1:number_cells_all
    neighbourhood.neurons(i_cell).params(batch_ID)=parameter_history(end,i_cell);
    neighbourhood.neurons(i_cell).posterior_stat(batch_ID)=...
        calculate_posterior(parameter_history(end,i_cell),quantile_prob);
end

%%
% %%
% traces=[];traces2=[];
% for i = 1:length(parameter_history)
% traces(i)=parameter_history(i).PR.mean;
% traces2(i)=parameter_history(i).PR.log_sigma;
% end
% figure(1)
% plot(traces)
% 
% figure(2)
% plot(traces2)
% figure(3)
% plot(elbo_rec(3:end))
% 
% %% Update the parameters in neighbourhood 
% 
% batch_ID=neighbourhood.batch_ID;
% quantile_prob=group_profile.regroup_func_params.quantile_prob;
% quantile_prob=[0.1 0.5 0.9];
% for i_batch = 1:batch_ID
%     for i_cell = 1:number_cells_all
%         if i_batch == 1
%         neighbourhood.neurons(i_cell).params(i_batch)=parameter_history(1,i_cell);
%         neighbourhood.neurons(i_cell).posterior_stat(i_batch)=...
%             calculate_posterior(parameter_history(1,i_cell),quantile_prob);
%         else
%             
%         neighbourhood.neurons(i_cell).params(i_batch)=parameter_history(end,i_cell);
%         neighbourhood.neurons(i_cell).posterior_stat(i_batch)=...
%             calculate_posterior(parameter_history(end,i_cell),quantile_prob);
%         end
%     end
% end
% %%
% neighbourhood.neurons(i_cell).posterior_stat(end).PR.mean
% neighbourhood.neurons(i_cell).posterior_stat(end).gain.mean
% neighbourhood.neurons(i_cell).posterior_stat(end).delay_mu.mean
% neighbourhood.neurons(i_cell).posterior_stat(end).delay_sigma.mean
% neurons.truth
% %%
% neighbourhood.neurons(i_cell).posterior_stat(end).PR.sd-neighbourhood.neurons(i_cell).posterior_stat(1).PR.sd
% neighbourhood.neurons(i_cell).posterior_stat(end).gain.sd-neighbourhood.neurons(i_cell).posterior_stat(1).gain.sd
% neighbourhood.neurons(i_cell).posterior_stat(end).delay_mu.sd-neighbourhood.neurons(i_cell).posterior_stat(1).delay_mu.sd
% neighbourhood.neurons(i_cell).posterior_stat(end).delay_sigma.sd-neighbourhood.neurons(i_cell).posterior_stat(1).delay_sigma.sd
% 
% %%
% simulation_params=experiment_setup.sim;
% est_shapes = neighbourhood.neurons(i_cell).posterior_stat(end).shapes;
% true_shape=neighbourhood.neurons(i_cell).truth.shape;
% this_locs=est_shapes.locations;
% this_size = griddata(simulation_params.mesh_grid(:,1),simulation_params.mesh_grid(:,2),simulation_params.mesh_grid(:,3),...
%     true_shape,this_locs(:,1),this_locs(:,2),this_locs(:,3));
% %%
% this_size*neighbourhood.neurons.truth.optical_gain
% est_shapes.mean*neighbourhood.neurons(i_cell).posterior_stat(end).gain.mean