function [logdist]=get_logdistribution(variational_samples,raw_samples, params)
%         logprior_shifts = sum(log( normpdf(shifts,0,sigma_shift)));
%         logprior_kernels=log( normpdf(log(taus),0, 2))+ log( normpdf(log(sigmasGP),0, 2));
%         logprior_sigma=log( normpdf(log(sigmas),0,2));
%%
epsilon=1e-4; % to prevent log(0)
fldnames = fieldnames(params(1));
n_cell = length(params);

logdist_mat = zeros(n_cell, length(fldnames));
for i_cell = 1:n_cell
    for i_field = 1:length(fldnames)
        this_sample=variational_samples(i_cell).(fldnames{i_field});
        this_raw_sample=raw_samples(i_cell).(fldnames{i_field});
        this_params=params(i_cell).(fldnames{i_field});
        switch this_params.dist
            case 'normal'
                this_logdist=log(normpdf(this_raw_sample,this_params.mean,exp(this_params.log_sigma)));
            case 'log-normal'
                this_logdist=log(normpdf(this_raw_sample,this_params.mean,exp(this_params.log_sigma)))+...
                    log(1/this_sample);
            case 'logit-normal'
                this_logdist=   log(normpdf(this_raw_sample,this_params.mean,exp(this_params.log_sigma))/...
                        ((this_sample-this_params.bounds.low)*(this_params.bounds.up-this_sample)) *(this_params.bounds.up-this_params.bounds.low));
                    
            case 'spiked-logit-normal'
                zero_prob =exp(this_params.prob_logit)/(1+exp(this_params.prob_logit));
                
                if this_samples == 0
                    this_logdist=   log(max(epsilon,zero_prob));
                else
                    this_logdist=   log(max(epsilon,1-zero_prob))+...
                    log(normpdf(this_raw_sample,this_params.mean,exp(this_params.log_sigma))/...
                        ((this_sample-this_params.bounds.low)*(this_params.bounds.up-this_sample)) *(this_params.bounds.up-this_params.bounds.low));
                    
                
                end
                             
        end
        logdist_mat(i_cell,i_field)=this_logdist;
        
        if (strcmp(this_params.type, 'common') & (i_cell >1)) %only count the common parameter once
            logdist_mat(i_cell,i_field)=0;
        end
        
    end
end
logdist=sum(sum(logdist_mat));
