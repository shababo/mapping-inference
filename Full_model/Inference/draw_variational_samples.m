function [parameter_sample, raw_sample] = draw_variational_samples(parameter_sample,...
    parameter_current,bounds)

%%
bounded_params_list= fieldnames(bounds);
n_cell=length(parameter_current);
params_list=fieldnames(parameter_current(1));
raw_sample=parameter_sample;
% maybe check if the parameters match between parameter_sample and
% parameter_current
S=size(parameter_sample.(params_list{1}),2);
for i_param = 1:length(params_list)
    for i_cell = 1:n_cell
    %%
    this_param=parameter_current(i_cell).(params_list{i_param});
    switch this_param.family
        case 'log-normal'
            [this_param_sample, this_raw_sample]=draw_lognormal(this_param, S);
        case 'logit-normal'
            [this_param_sample, this_raw_sample]=draw_logitnormal(this_param, S);
            if ismember(params_list{i_param}, bounded_params_list)
                this_param_sample= this_param_sample*range(bounds.(params_list{i_param}))+bounds.(params_list{i_param})(1);
            end
        case 'spiked-logit-normal'
            % placeholder
    end
    % if there are bounds to be enforced 
    
    parameter_sample.(params_list{i_param})(i_cell,:)=this_param_sample;
    raw_sample.(params_list{i_param})(i_cell,:)=this_raw_sample;
    
    end
end




