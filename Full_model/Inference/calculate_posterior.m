function [posterior_params] = calculate_posterior(parameters,gamma_bound,gain_bound,quantiles_prob,varargin)


if ~isempty(varargin) && ~isempty(varargin{1})
  spike_indicator=varargin{1};
else
  spike_indicator=false;
end

n_cell=length(parameters);
% Initialize storage for the fitted parameters in the experiment
[mean_gamma_temp, ~] = calculate_posterior_mean(...
    [parameters(:).alpha]',exp([parameters(:).beta]'),gamma_bound.low,gamma_bound.up);
if spike_indicator
    v_pi=exp([parameters(:).p_logit]')./(1+exp([parameters(:).p_logit]'));
    mean_gamma_temp=mean_gamma_temp.*(1-v_pi);
else
    v_pi=zeros(n_cell,1);
end
[mean_gain_temp, var_gain_temp] = calculate_posterior_mean(...
    [parameters(:).alpha_gain]',exp([parameters(:).beta_gain]'),gain_bound.low,gain_bound.up);
[low_gamma_temp, up_gamma_temp]=calculate_posterior_quatiles(quantiles_prob,...
    [parameters(:).alpha]',exp([parameters(:).beta]'),gamma_bound.low,gamma_bound.up);


posterior_params=struct([]);
temp=num2cell(mean_gamma_temp);[posterior_params(1:n_cell).gamma_mean]=temp{:};
temp=num2cell(low_gamma_temp);[posterior_params(1:n_cell).gamma_lower_quantile]=temp{:};
temp=num2cell(up_gamma_temp);[posterior_params(1:n_cell).gamma_upper_quantile]=temp{:};
temp=num2cell(1-v_pi);[posterior_params(1:n_cell).nonzero_prob]=temp{:};
temp=num2cell(mean_gain_temp);[posterior_params(1:n_cell).gain_mean]=temp{:};
temp=num2cell(var_gain_temp);[posterior_params(1:n_cell).gain_var]=temp{:};



      
        
end




