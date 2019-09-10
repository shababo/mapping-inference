function [params,change]=incorporate_gradient(params, new_gradient,varargin)
%%

if ~isempty(varargin)
    marginal_flag = varargin{1};
else
    marginal_flag = false;
end
fldnames = fieldnames(new_gradient(1,1));
n_cell = length(params);
for i_cell = 1:n_cell
    for i_field = 1:length(fldnames)
%         if marginal_flag & ( strcmp(fldnames{i_field},'shapes') | strcmp(fldnames{i_field},'delay_mu') | strcmp(fldnames{i_field},'delay_sigma'))
%         else
        dims=size(params(i_cell).(fldnames{i_field}).mean);
        params(i_cell).(fldnames{i_field}).mean=...
            params(i_cell).(fldnames{i_field}).mean+reshape(new_gradient(i_cell).(fldnames{i_field}).mean,dims);
        params(i_cell).(fldnames{i_field}).log_sigma=...
            params(i_cell).(fldnames{i_field}).log_sigma+reshape(new_gradient(i_cell).(fldnames{i_field}).sigma,dims);
        if strcmp(params(i_cell).(fldnames{i_field}).mean,  'spiked-logit-normal')
              params(i_cell).(fldnames{i_field}).prob_logit=...
            params(i_cell).(fldnames{i_field}).prob_logit+new_gradient(i_cell).(fldnames{i_field}).prob_logit;
        end
%         end
        if ( strcmp(fldnames{i_field},'shapes') | strcmp(fldnames{i_field},'xy') | strcmp(fldnames{i_field},'z'))  & strcmp( params(i_cell).(fldnames{i_field}).dist,'mvn')
            
       %Dmat= diag(diag(params(i_cell).(fldnames{i_field}).Sigma_inv).*exp(-params(i_cell).(fldnames{i_field}).log_sigma) );
       Dmat= diag(exp(-params(i_cell).(fldnames{i_field}).log_sigma) );
       params(i_cell).(fldnames{i_field}).Sigma_tilde_inv=params(i_cell).(fldnames{i_field}).Sigma_inv+Dmat;
        params(i_cell).(fldnames{i_field}).Sigma_tilde=inv(params(i_cell).(fldnames{i_field}).Sigma_tilde_inv);
        params(i_cell).(fldnames{i_field}).Sigma_tilde=(params(i_cell).(fldnames{i_field}).Sigma_tilde+params(i_cell).(fldnames{i_field}).Sigma_tilde')/2;
  
            
            
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
        mag_param_mean=mag_param_mean+sum(abs(params(i_cell).(fldnames{i_field}).mean));
        mag_param_sigma=mag_param_sigma+sum(abs(params(i_cell).(fldnames{i_field}).log_sigma));
        
        mag_grad_mean=mag_grad_mean+sum(abs(new_gradient(i_cell).(fldnames{i_field}).mean));
        mag_grad_sigma=mag_grad_sigma+sum(abs(new_gradient(i_cell).(fldnames{i_field}).sigma));
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
