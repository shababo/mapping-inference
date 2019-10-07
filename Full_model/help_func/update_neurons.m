function [neurons]=update_neurons(new_ID,neurons,est_parameter,inference_params)
%% Updating the fields of neurons using the estimated parameters 
for i_cell = 1:length(neurons)
    neurons(i_cell).params(new_ID)=est_parameter(i_cell);
    neurons(i_cell).posterior_stat(new_ID)=...
        calculate_posterior(est_parameter(i_cell),inference_params.quantile_prob);
%     neurons(i_cell).prior_stat(new_ID)=calculate_posterior(parameter_history(1,i_cell),inference_params.quantile_prob);
end

