function [this_gradient]=get_variational_gradient(variational_samples,raw_samples,params)

% calculate the gradients
fldnames = fieldnames(params(1));
n_cell = length(params);
this_gradient(n_cell)= struct;
for i_cell = 1:n_cell
    for i_field = 1:length(fldnames)
        this_params=params(i_cell).(fldnames{i_field});
        this_sample=variational_samples(i_cell).(fldnames{i_field});
        this_raw_sample = raw_samples(i_cell).(fldnames{i_field});
        if ~strcmp(fldnames{i_field},'shapes')
        switch this_params.dist
            case {'normal','log-normal', 'logit-normal'}
                dmean =   -(this_params.mean-this_raw_sample)./exp(2*this_params.log_sigma);
                dsigma = -1+ (this_params.mean-this_raw_sample).^2./exp(2*this_params.log_sigma);
            case {'spiked-logit-normal'}
                dlogit= (this_sample==0)/(1+exp(this_params.prob_logit))-...
                    (this_sample==0)*exp(this_params.prob_logit)/(1+exp(this_params.prob_logit));
                dmean =   -(this_params.mean-this_raw_sample)./exp(2*this_params.log_sigma);
                dsigma = -1+ (this_params.mean-this_raw_sample).^2./exp(2*this_params.log_sigma);
            this_gradient(i_cell).(fldnames{i_field}).prob_logit=dlogit;
        
        end
        else % shapes
            for i_loc = 1:length(this_raw_sample)
                switch this_params.dist
                    case {'normal','log-normal', 'logit-normal'}
                        dmean(i_loc) =   -(this_params.mean(i_loc)-this_raw_sample(i_loc))./exp(2*this_params.log_sigma(i_loc));
                        dsigma(i_loc) = -1+ (this_params.mean(i_loc)-this_raw_sample(i_loc)).^2./exp(2*this_params.log_sigma(i_loc));    
                end
            end
        end
        
        this_gradient(i_cell).(fldnames{i_field}).type =  this_params.type;
        this_gradient(i_cell).(fldnames{i_field}).mean = dmean;
        this_gradient(i_cell).(fldnames{i_field}).sigma = dsigma;
        
    end
end


