function [this_gradient]=get_variational_gradient(variational_sample,params)

% calculate the gradients
fldnames = fieldnames(params(1));
n_cell = length(params);
this_gradient(n_cell)= struct;
for i_cell = 1:n_cell
    for i_field = 1:length(fldnames)
        this_params=params(i_cell).(fldnames{i_field});
        this_sample=variational_sample(i_cell).(fldnames{i_field});
        
        if strcmp(this_params.dist,'log-normal')
           this_sample = log(this_sample); 
        end
        
        switch this_params.dist
            case {'normal','log-normal'}
            dmean =   -(this_params.mean-this_sample)./exp(2*this_params.log_sigma);
            dsigma = -1+ (this_params.mean-this_sample).^2./exp(2*this_params.log_sigma);
        end
        this_gradient(i_cell).(fldnames{i_field}).type =  this_params.type;
        
        this_gradient(i_cell).(fldnames{i_field}).mean = dmean;
        this_gradient(i_cell).(fldnames{i_field}).sigma = dsigma;
        
    end
end


