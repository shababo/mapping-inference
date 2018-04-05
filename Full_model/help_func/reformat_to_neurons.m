function [output]=reformat_to_neurons(fitted_params,property_type,dist_type)

% fitted_params=parameter_history(end,1);
% property_type='gamma';
output=struct;
switch property_type
    case 'PR'
        
        output.type=dist_type;
        output.pi_logit=fitted_params.p_logit;
        output.alpha=fitted_params.alpha;
        output.beta=fitted_params.beta;
    case 'gain'
        output.type=dist_type;
        output.pi_logit=-Inf;
        output.alpha=fitted_params.alpha_gain;
        output.beta=fitted_params.beta_gain;
    case 'delay_mu'
        output.type=dist_type;
        output.pi_logit=-Inf;
        output.alpha=fitted_params.alpha_m;
        output.beta=fitted_params.beta_m;
    case 'delay_sigma'
        output.type=dist_type;
        output.pi_logit=-Inf;
        output.alpha=fitted_params.alpha_s;
        output.beta=fitted_params.beta_s;
end
end
