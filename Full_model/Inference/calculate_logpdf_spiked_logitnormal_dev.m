function [logpdf] = calculate_logpdf_spiked_logitnormal_dev(...
    parameter_dist,logit_gamma,gamma_sample,logit_gain,gain_sample,...
    logit_delay_mu,delay_mu_sample,logit_delay_sigma,delay_sigma_sample,...
    delay_mu_bound,delay_sigma_bound, gamma_bound,gain_bound,varargin)



if ~isempty(varargin) && ~isempty(varargin{1})
    spike_indicator=varargin{1};
else
    spike_indicator=false;
end

if spike_indicator
    v_pi=exp([parameter_dist(:).p_logit]')./(1+exp([parameter_dist(:).p_logit]'));
    

    logpdf=    log(max(0.001,v_pi)).*(gamma_sample==0)+...
        log(max(0.001, 1-v_pi)).*(gamma_sample>0).*log( min(1000,max(0.0001,...
        normpdf(logit_gamma,[parameter_dist(:).alpha]',exp([parameter_dist(:).beta]'))...
        ./( (gamma_sample-gamma_bound.low ).*(gamma_bound.up-gamma_sample)) *(gamma_bound.up-gamma_bound.low))))...
        +(gamma_sample>0).*log( min(1000,max(0.0001,...
        normpdf(logit_gain,[parameter_dist(:).alpha_gain]',exp([parameter_dist(:).beta_gain]'))...
        ./( (gain_sample-gain_bound.low ).*(gain_bound.up-gain_sample)) *(gain_bound.up-gain_bound.low))))...
        + log( min(1000,max(0.0001,...
        normpdf(logit_delay_mu,[parameter_dist(:).alpha_m]',exp([parameter_dist(:).beta_m]'))...
        ./( (delay_mu_sample-delay_mu_bound.low ).*(delay_mu_bound.up-delay_mu_sample)) *(delay_mu_bound.up-delay_mu_bound.low))))...
        + log( min(1000,max(0.0001,...
        normpdf(logit_delay_sigma,[parameter_dist(:).alpha_s]',exp([parameter_dist(:).beta_s]'))...
        ./( (delay_sigma_sample-delay_sigma_bound.low ).*(delay_sigma_bound.up-delay_sigma_sample)) *(delay_sigma_bound.up-delay_sigma_bound.low))));
else

    logpdf=   log( min(1000,max(0.0001,...
        normpdf(logit_gamma,[parameter_dist(:).alpha]',exp([parameter_dist(:).beta]'))...
        ./( (gamma_sample-gamma_bound.low ).*(gamma_bound.up-gamma_sample)) *(gamma_bound.up-gamma_bound.low))))...
        +log( min(1000,max(0.0001,...
        normpdf(logit_gain,[parameter_dist(:).alpha_gain]',exp([parameter_dist(:).beta_gain]'))...
        ./( (gain_sample-gain_bound.low ).*(gain_bound.up-gain_sample)) *(gain_bound.up-gain_bound.low))))...
        + log( min(1000,max(0.0001,...
        normpdf(logit_delay_mu,[parameter_dist(:).alpha_m]',exp([parameter_dist(:).beta_m]'))...
        ./( (delay_mu_sample-delay_mu_bound.low ).*(delay_mu_bound.up-delay_mu_sample)) *(delay_mu_bound.up-delay_mu_bound.low))))...
        + log( min(1000,max(0.0001,...
        normpdf(logit_delay_sigma,[parameter_dist(:).alpha_s]',exp([parameter_dist(:).beta_s]'))...
        ./( (delay_sigma_sample-delay_sigma_bound.low ).*(delay_sigma_bound.up-delay_sigma_sample)) *(delay_sigma_bound.up-delay_sigma_bound.low))));
    
end

end