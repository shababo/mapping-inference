function [params,change]=incorporate_gradient(params, new_gradient)

fldnames = fieldnames(new_gradient(1,1));
n_cell = length(params);


for i_cell = 1:n_cell
    for i_field = 1:length(fldnames)
        params(i_cell).(fldnames{i_field}).mean=...
            params(i_cell).(fldnames{i_field}).mean+new_gradient(i_cell).(fldnames{i_field}).mean;
        params(i_cell).(fldnames{i_field}).log_sigma=...
            params(i_cell).(fldnames{i_field}).log_sigma+new_gradient(i_cell).(fldnames{i_field}).sigma;
        if strcmp(params(i_cell).(fldnames{i_field}).mean,  'spiked-logit-normal')
              params(i_cell).(fldnames{i_field}).prob_logit=...
            params(i_cell).(fldnames{i_field}).prob_logit+new_gradient(i_cell).(fldnames{i_field}).prob_logit;
      
        end
    end
end



all_changes=zeros(length(fldnames),1);
for i_field = 1:length(fldnames)
    mag_param_mean = 0; mag_grad_mean = 0;
    mag_param_sigma = 0; mag_grad_sigma = 0;
     if strcmp(params(i_cell).(fldnames{i_field}).mean,  'spiked-logit-normal')
        mag_param_prob_logit = 0; mag_grad_prob_logit = 0;
   
     end
    for i_cell = 1:n_cell
        mag_param_mean=mag_param_mean+abs(params(i_cell).(fldnames{i_field}).mean);
        mag_param_sigma=mag_param_sigma+abs(params(i_cell).(fldnames{i_field}).log_sigma);
        
        mag_grad_mean=mag_grad_mean+abs(new_gradient(i_cell).(fldnames{i_field}).mean);
        mag_grad_sigma=mag_grad_sigma+abs(new_gradient(i_cell).(fldnames{i_field}).sigma);
         if strcmp(params(i_cell).(fldnames{i_field}).mean,  'spiked-logit-normal')
         mag_param_prob_logit=mag_param_prob_logit+abs(params(i_cell).(fldnames{i_field}).prob_logit);
        mag_grad_prob_logit=mag_grad_prob_logit+abs(new_gradient(i_cell).(fldnames{i_field}).prob_logit);
       
         end
    end
    all_changes(i_field)= mag_grad_mean/mag_param_mean+ mag_grad_sigma/mag_param_sigma;
     if strcmp(params(i_cell).(fldnames{i_field}).mean,  'spiked-logit-normal')
     all_changes(i_field)= all_changes(i_field)+ mag_grad_prob_logit/mag_param_prob_logit;
     end
end

 change= sum(all_changes);
