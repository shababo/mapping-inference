function [variational_samples,raw_samples] = draw_samples_from_var_dist(variational_params)
%%
fldnames = fieldnames(variational_params(1));


n_cell = length(variational_params);
clear('variational_samples')
variational_samples(n_cell)=struct;
raw_samples(n_cell)=struct;
for i_cell = 1:n_cell
    for i_field = 1:length(fldnames)
        this_params=variational_params(i_cell).(fldnames{i_field});
%         this_sample=struct;
        switch this_params.dist
            case 'normal'
                
                this_raw_sample=normrnd(this_params.mean,exp(this_params.log_sigma));
                this_sample=this_raw_sample;
            case 'log-normal'
                this_raw_sample=normrnd(this_params.mean,exp(this_params.log_sigma));
                this_sample=exp(this_raw_sample);
                
            case 'logit-normal'
               this_raw_sample=normrnd(this_params.mean,exp(this_params.log_sigma));
               this_sample = exp(this_raw_sample)./(1+exp(this_raw_sample))*...
                   (this_params.bounds.up-this_params.bounds.low)+this_params.bounds.low;
                
            case 'spiked-logit-normal'
           this_raw_sample=normrnd(this_params.mean,exp(this_params.log_sigma));
               this_sample = exp(this_raw_sample)./(1+exp(this_raw_sample))*...
                   (this_params.bounds.up-this_params.bounds.low)+this_params.bounds.low;
               zero_prob = exp(this_params.prob_logit)./(exp(this_params.prob_logit)+1);
               this_spike =   rand(1) > (zero_prob);
                
               this_sample = this_spike*this_sample;
        end
        
        variational_samples(i_cell).(fldnames{i_field}) = this_sample;
        raw_samples(i_cell).(fldnames{i_field}) = this_raw_sample;
        
        if strcmp(this_params.type, 'common')
            variational_samples(i_cell).(fldnames{i_field})  =  variational_samples(1).(fldnames{i_field});
            raw_samples(i_cell).(fldnames{i_field})  =  raw_samples(1).(fldnames{i_field});
        
        end
        
    end
end
