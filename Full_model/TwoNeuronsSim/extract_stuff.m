function [this_params_mean]=extract_stuff(current_parameter,all_bounds,quantile_prob,...
        this_field)
    current_params=reformat_to_neurons(current_parameter,this_field,'spiked_logit_normal');
    %     group_profile=experiment_setup.groups.(this_neighbourhood.neurons(i_cell).group_ID);
    bounds= all_bounds.(this_field);
   this_params=...
        calculate_posterior(current_params,bounds,quantile_prob);
    this_params_mean=this_params.mean;
    