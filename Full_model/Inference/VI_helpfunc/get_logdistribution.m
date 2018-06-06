function [logdist]=get_logdistribution(variational_sample,params)
%         logprior_shifts = sum(log( normpdf(shifts,0,sigma_shift)));
%         logprior_kernels=log( normpdf(log(taus),0, 2))+ log( normpdf(log(sigmasGP),0, 2));
%         logprior_sigma=log( normpdf(log(sigmas),0,2));
     %%    
fldnames = fieldnames(params(1));
n_cell = length(params);

logdist_mat = zeros(n_cell, length(fldnames));
for i_cell = 1:n_cell
    for i_field = 1:length(fldnames)
        this_sample=variational_sample(i_cell).(fldnames{i_field});
        
        this_params=params(i_cell).(fldnames{i_field});
        switch this_params.dist
            case 'normal'
                this_logdist=log(normpdf(this_sample,this_params.mean,exp(this_params.log_sigma)));
            case 'log-normal'
                   this_logdist=log(normpdf(log(this_sample),this_params.mean,exp(this_params.log_sigma)));
        end
        logdist_mat(i_cell,i_field)=this_logdist; 
        
        if (strcmp(this_params.type, 'common') & (i_cell >1)) %only count the common parameter once
            logdist_mat(i_cell,i_field)=0;
        end
        
    end
end
logdist=sum(sum(logdist_mat));
