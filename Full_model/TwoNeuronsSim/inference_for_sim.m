function [fitted_neurons]=inference_for_sim(mpp,neurons,bg_rate,z_dictionary,shape_var)

%% Preparation
experiment_setup=get_experiment_setup('szchen-sim');
group_profile=experiment_setup.groups.undefined;
inference_params=group_profile.inference_params;
% inference_params.likelihood=@calculate_likelihood_intensity_for_VI;
inference_params.likelihood=@lif_glm_firstspike_loglikelihood_for_VI;
inference_params.bounds.PR=[0.01 1];
inference_params.bounds.gain=[0.005 0.06];
% inference_params.MCsamples_for_gradient=100;
inference_params.convergence_threshold=1e-4;
inference_params.maxit=2000;
% inference_params.step_size_max=1;

prior_info=experiment_setup.prior_info;

%% Initialization:
temp_output=struct;
temp_output.PR_params=struct;
temp_output.PR_params.pi_logit=-Inf*ones(2,1);
temp_output.PR_params.alpha=0*ones(2,1);
temp_output.PR_params.beta=1*ones(2,1);

temp_output.gain_params=temp_output.PR_params;
temp_output.delay_mu_params=temp_output.PR_params;
temp_output.delay_sigma_params=temp_output.PR_params;

%% Formatting
number_of_stim_cells=length(neurons);
variational_params=struct;
temp=num2cell(temp_output.PR_params.pi_logit); [variational_params(1:number_of_stim_cells).p_logit]=temp{:};
temp=num2cell(temp_output.PR_params.alpha); [variational_params(1:number_of_stim_cells).alpha]=temp{:};
temp=num2cell(temp_output.PR_params.beta); [variational_params(1:number_of_stim_cells).beta]=temp{:};
temp=num2cell(temp_output.gain_params.alpha); [variational_params(1:number_of_stim_cells).alpha_gain]=temp{:};
temp=num2cell(temp_output.gain_params.beta); [variational_params(1:number_of_stim_cells).beta_gain]=temp{:};
temp=num2cell(temp_output.delay_mu_params.alpha); [variational_params(1:number_of_stim_cells).alpha_m]=temp{:};
temp=num2cell(temp_output.delay_mu_params.beta); [variational_params(1:number_of_stim_cells).beta_m]=temp{:};
temp=num2cell(temp_output.delay_sigma_params.alpha); [variational_params(1:number_of_stim_cells).alpha_s]=temp{:};
temp=num2cell(temp_output.delay_sigma_params.beta); [variational_params(1:number_of_stim_cells).beta_s]=temp{:};

prior_params=variational_params;

%% SVI:
switch shape_var
    case 'mean_shape'
        stim_size=reshape([mpp(:).stimulation], [2 length(mpp)])';
        [parameter_history,lklh_rec] = fit_VI(...
            stim_size, mpp, bg_rate,...
            variational_params,prior_params,...
            inference_params,prior_info);
    case 'dic_shape'
        [parameter_history,lklh_rec] = fit_VI_dev(...
             mpp, bg_rate,z_dictionary,...
            variational_params,prior_params,...
            inference_params,prior_info);
end

%% Reformat the output 
clear('fitted_neurons')
quantile_prob=group_profile.regroup_func_params.quantile_prob;
 stim_size=reshape([mpp(:).stimulation], [2 length(mpp)])';
       
fitted_neurons(length(neurons))=struct;
for this_cell = 1:length(neurons)
    fitted_neurons(this_cell).delay_var=extract_stuff(parameter_history(end, this_cell),group_profile.inference_params.bounds,quantile_prob,...
        'delay_sigma')^2;
    fitted_neurons(this_cell).delay_mean=extract_stuff(parameter_history(end, this_cell),group_profile.inference_params.bounds,quantile_prob,...
        'delay_mu');
    fitted_neurons(this_cell).gain=extract_stuff(parameter_history(end, this_cell),group_profile.inference_params.bounds,quantile_prob,...
        'gain');
    fitted_neurons(this_cell).PR=extract_stuff(parameter_history(end, this_cell),group_profile.inference_params.bounds,quantile_prob,...
        'PR');
    fitted_neurons(this_cell).mpp=mpp;
    fitted_neurons(this_cell).stim=stim_size(:,this_cell);
    fitted_neurons(this_cell).location=neurons(this_cell).location;
end
%% Visualize the path 
figure(1)
plot(lklh_rec)

%%
figure(2)
plot([parameter_history(:,2).alpha])