function [logdist]=get_logdistribution(variational_samples,raw_samples, params,varargin)
%         logprior_shifts = sum(log( normpdf(shifts,0,sigma_shift)));
%         logprior_kernels=log( normpdf(log(taus),0, 2))+ log( normpdf(log(sigmasGP),0, 2));
%         logprior_sigma=log( normpdf(log(sigmas),0,2));
%%

% if ~isempty(varargin)
%    only_PR=varargin{1};
% else
%    only_PR=false;
% end

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
                this_logdist=log(normpdf(this_raw_sample,this_params.mean,exp(this_params.log_sigma/2)));
            case 'log-normal'
                this_logdist=log(normpdf(this_raw_sample,this_params.mean,exp(this_params.log_sigma/2)))+...
                    log(1/this_sample);
            case 'logit-normal'
                this_logdist=   log(normpdf(this_raw_sample,this_params.mean,exp(this_params.log_sigma/2))/...
                    ((this_sample-this_params.bounds.low+epsilon)*(this_params.bounds.up-this_sample+epsilon)) *(this_params.bounds.up-this_params.bounds.low));
            case 'spiked-logit-normal'
                zero_prob =exp(this_params.prob_logit)/(1+exp(this_params.prob_logit/2));
                
                if this_samples == 0
                    this_logdist=   log(max(epsilon,zero_prob));
                else
                    this_logdist=   log(max(epsilon,1-zero_prob))+...
                        log(normpdf(this_raw_sample,this_params.mean,exp(this_params.log_sigma))/...
                        ((this_sample-this_params.bounds.low)*(this_params.bounds.up-this_sample)) *(this_params.bounds.up-this_params.bounds.low));
                end
            case 'mvn-logit'
                % transfer the mean:
                this_raw_sample=raw_samples(i_cell).(fldnames{i_field});
                this_params=params(i_cell).(fldnames{i_field});
                
                this_mean=(this_params.bounds.up-this_params.bounds.low).*exp(this_params.mean)./(1+exp(this_params.mean)) +this_params.bounds.low;
                %this_mean=this_params.mean;
                this_logdist=log(mvnpdf(reshape(this_raw_sample, size(this_mean)),this_mean,this_params.Sigma_tilde));
            
            case 'mvn'
                this_raw_sample=raw_samples(i_cell).(fldnames{i_field});
                this_params=params(i_cell).(fldnames{i_field});
                
                this_mean=this_params.mean;
                %this_mean=this_params.mean;
                this_logdist=log(mvnpdf(reshape(this_raw_sample, size(this_mean)),this_mean,this_params.Sigma_tilde));
        end
        
        if isnan(this_logdist) | isinf(this_logdist)
            this_logdist=log(epsilon);
        end
        logdist_mat(i_cell,i_field)=this_logdist;
        
        if (strcmp(this_params.type, 'common') & (i_cell >1)) %only count the common parameter once
            logdist_mat(i_cell,i_field)=0;
        end
        
    end
end
% if only_PR
%     i_field = 2;
% logdist=sum(sum(logdist_mat(:,i_field)));
%
% else
logdist=sum(sum(logdist_mat));
% end
