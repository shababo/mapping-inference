function [variational_sample] = draw_samples_from_var_dist(variational_params)
%%    
fldnames = fieldnames(variational_params(1));


n_cell = length(variational_params);
clear('variational_sample')
variational_sample(n_cell)=struct;
for i_cell = 1:n_cell
    for i_field = 1:length(fldnames)
        this_params=variational_params(i_cell).(fldnames{i_field});
        this_sample=struct;
        switch this_params.dist
            case 'normal'
                this_sample=normrnd(this_params.mean,exp(this_params.log_sigma));
            case 'log-normal'
                this_sample=exp(normrnd(this_params.mean,exp(this_params.log_sigma)));
        end
        
        variational_sample(i_cell).(fldnames{i_field}) = this_sample;
        
        if strcmp(this_params.type, 'common')
            variational_sample(i_cell).(fldnames{i_field})  =  variational_sample(1).(fldnames{i_field});
        end
        
    end
end
     