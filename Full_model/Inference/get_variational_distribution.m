function [] = get_variational_distribution(...
    parameter_distribution,parameter_sample,raw_sample)


% for each cell and each parameter from the MC samples
% calculate the variational distribution
% and calculate the derivatives w.r.t. 

n_cell=length(parameter_distribution);
params_list=fieldnames(parameter_distribution(1));
S=size(parameter_sample.(params_list{1}),2);
log_dist = struct;
% initial all the fields:
for i_param = 1:length(params_list)
    
end 
for i_cell = 1:n_cell    
    for s=1:S
        % Calculate and store the log p for all parameter
        % Calculate the gradient too 
        
        % need to define a new structure to save it.
        for i_param = 1:length(params_list)
             
        end
    end
end

if spike_indicator
    v_pi=exp([parameter_dist(:).p_logit]')./(1+exp([parameter_dist(:).p_logit]'));
    

    logpdf=    log(max(0.001,v_pi)).*(gamma_sample==0)+...
        log(max(0.001, 1-v_pi)).*(gamma_sample>0).*log( min(1000,max(0.0001,...
        normpdf(logit_gamma,[parameter_dist(:).alpha]',exp([parameter_dist(:).beta]'))...
        ./( (gamma_sample-gamma_bound.low ).*(gamma_bound.up-gamma_sample)) *(gamma_bound.up-gamma_bound.low))))...
        +(gamma_sample>0).*log( min(1000,max(0.0001,...
        normpdf(logit_gain,[parameter_dist(:).alpha_gain]',exp([parameter_dist(:).beta_gain]'))...
        ./( (gain_sample-gain_bound.low ).*(gain_bound.up-gain_sample)) *(gain_bound.up-gain_bound.low))));
    
else

    logpdf=   log( min(1000,max(0.0001,...
        normpdf(logit_gamma,[parameter_dist(:).alpha]',exp([parameter_dist(:).beta]'))...
        ./( (gamma_sample-gamma_bound.low ).*(gamma_bound.up-gamma_sample)) *(gamma_bound.up-gamma_bound.low))))...
        +log( min(1000,max(0.0001,...
        normpdf(logit_gain,[parameter_dist(:).alpha_gain]',exp([parameter_dist(:).beta_gain]'))...
        ./( (gain_sample-gain_bound.low ).*(gain_bound.up-gain_sample)) *(gain_bound.up-gain_bound.low))));
    
    
end