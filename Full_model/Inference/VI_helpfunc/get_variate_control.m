function [this_gradient]=get_variate_control(lklhweights,this_gradient)

fldnames = fieldnames(this_gradient(1));
n_cell = length(this_gradient);
for i_cell = 1:n_cell
    for i_field = 1:length(fldnames)
        
        this_gradient(i_cell).(fldnames{i_field}).mean_f=...
            this_gradient(i_cell).(fldnames{i_field}).mean*lklhweights;
        this_gradient(i_cell).(fldnames{i_field}).mean_h=...
            this_gradient(i_cell).(fldnames{i_field}).mean;
        this_gradient(i_cell).(fldnames{i_field}).sigma_f=...
            this_gradient(i_cell).(fldnames{i_field}).sigma*lklhweights;
        this_gradient(i_cell).(fldnames{i_field}).sigma_h=...
            this_gradient(i_cell).(fldnames{i_field}).sigma;
        if strcmp(this_gradient(i_cell).(fldnames{i_field}).type, 'spiked-logit-normal')
            this_gradient(i_cell).(fldnames{i_field}).prob_logit_f=...
                this_gradient(i_cell).(fldnames{i_field}).prob_logit*lklhweights;
        end
    end
end


