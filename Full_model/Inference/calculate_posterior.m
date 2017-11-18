function [current_params] = calculate_posterior(current_params,bounds,quantile_prob)
              
switch current_params.type
    case 'spiked_logit_normal'
        

n_cell=length(current_params);
% Initialize storage for the fitted parameters in the experiment
[mean_temp, variance_temp] = calculate_posterior_mean(...
    current_params.alpha,exp(current_params.beta),bounds(1),bounds(2));


nonzero_prob_temp=1-exp( current_params.pi_logit)./(1+exp(current_params.pi_logit));
mean_temp=mean_temp.*nonzero_prob_temp;
[low_temp, up_temp]=calculate_posterior_quatiles(quantile_prob,...
    current_params.alpha,exp(current_params.beta),bounds(1),bounds(2));

current_params.mean=mean_temp;
current_params.variance=variance_temp;
current_params.upper_quantile=up_temp;
current_params.lower_quantile=low_temp;
current_params.nonzero_prob=nonzero_prob_temp;


% posterior_params=struct([]);
% temp=num2cell(mean_gamma_temp);[posterior_params(1:n_cell).gamma_mean]=temp{:};
% temp=num2cell(low_gamma_temp);[posterior_params(1:n_cell).gamma_lower_quantile]=temp{:};
% temp=num2cell(up_gamma_temp);[posterior_params(1:n_cell).gamma_upper_quantile]=temp{:};
% temp=num2cell(1-v_pi);[posterior_params(1:n_cell).nonzero_prob]=temp{:};
% temp=num2cell(mean_gain_temp);[posterior_params(1:n_cell).gain_mean]=temp{:};
% temp=num2cell(var_gain_temp);[posterior_params(1:n_cell).gain_var]=temp{:};


    otherwise
        
end

      
        
end




