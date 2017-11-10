function [this_neighbourhood]=inference_undefined(...
    this_experiment_query,this_neighbourhood,experiment_setup)


        indicators_remained = find(ismember([mpp_undefined(:).batch],iter-(0:num_trace_back) ));
        mpp_remained=mpp_undefined(indicators_remained);
        trials_locations=reshape([mpp_remained(:).locations],n_spots_per_trial,[])';
        trials_powers=reshape([mpp_remained(:).power],n_spots_per_trial,[])';
        designs_remained = get_stim_size(pi_target,trials_locations,trials_powers);
           
        % include all cells that have been stimulated:
        cell_list= find(sum(designs_remained>stim_threshold,1)>0);
        designs_remained=designs_remained(:,cell_list);
        % Update variational and prior distribution
        variational_params=variational_params_path(iter,cell_list);
        prior_params=variational_params_path(max(iter-num_trace_back,1),cell_list);
        switch model_type
            case 1% working model
                lklh_func=@calculate_loglikelihood_bernoulli;
            case 2% working model with spike & slab
                lklh_func=@calculate_loglikelihood_bernoulli;
            case 3 % full model with first spike
                lklh_func=@lif_glm_firstspike_loglikelihood_for_VI;
            case 4 % full model with first event 
                lklh_func=@lif_glm_firstevent_loglikelihood_for_VI;
            otherwise
        end
        [parameter_history] = fit_VI(...
            designs_remained, mpp_remained, background_rate, ...
            prob_trace_full,stim_scale,minimum_stim_threshold,...
            variational_params,prior_params,gamma_bound,gain_bound,...
            S,epsilon,eta,eta_max,maxit,lklh_func,spike_indicator);
        
        variational_params_path(iter+1,cell_list)=parameter_history(end,:);
        
        [parameter_temp] = calculate_posterior(parameter_history(end,:),gamma_bound,gain_bound,quantiles_prob,spike_indicator);
        parameter_path(iter+1,cell_list)=parameter_temp;
        
end
