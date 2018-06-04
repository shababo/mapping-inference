function [new_gradient]=sum_gradient(gradients,eta,eta_max,iteration)

step_size =  eta*(iteration)^(-1/1.5);
%     step_size =  eta*(iteration)^(-1);
if eta_max>step_size
    eta_max=step_size;
end

fldnames = fieldnames(gradients(1,1));
n_cell = size(gradients,2);
S=size(gradients,1);
new_gradient=struct;
grad_max=0;
for i_cell = 1:n_cell
    for i_field = 1:length(fldnames)
        this_field_sigma_f=zeros(S,1);
        this_field_sigma_h=zeros(S,1);
        this_field_mean_f=zeros(S,1);
        this_field_mean_h=zeros(S,1);
        
        for s = 1:S
            this_field_mean_f(s) = gradients(s,i_cell).(fldnames{i_field}).mean_f;
            this_field_mean_h(s) = gradients(s,i_cell).(fldnames{i_field}).mean_h;
            
            this_field_sigma_f(s) = gradients(s,i_cell).(fldnames{i_field}).sigma_f;
            this_field_sigma_h(s) = gradients(s,i_cell).(fldnames{i_field}).sigma_h;
        end
        new_gradient(i_cell).(fldnames{i_field})=struct;
        new_gradient(i_cell).(fldnames{i_field}).type=gradients(1,i_cell).(fldnames{i_field}).type;
        
        if strcmp(gradients(1,i_cell).(fldnames{i_field}).type,'individual')
         a= quick_cov( this_field_mean_f, this_field_mean_h)+quick_cov( this_field_sigma_f, this_field_sigma_h);
         a=a/(quick_cov( this_field_mean_h, this_field_mean_h)+quick_cov( this_field_sigma_h, this_field_sigma_h));
         new_gradient(i_cell).(fldnames{i_field}).mean= ...
             step_size*(mean(this_field_mean_f)-a*mean(this_field_mean_h));
         new_gradient(i_cell).(fldnames{i_field}).sigma= ...
             step_size*(mean(this_field_sigma_f)-a*mean(this_field_sigma_h));
         
        else % if it is a common parameter: 
            new_gradient(i_cell).(fldnames{i_field}).mean= ...
             step_size*mean(this_field_mean_f);
         new_gradient(i_cell).(fldnames{i_field}).sigma= ...
             step_size*mean(this_field_sigma_f);
         
        end
        grad_max=max([grad_max abs( new_gradient(i_cell).(fldnames{i_field}).mean) abs( new_gradient(i_cell).(fldnames{i_field}).sigma)] );
    end
end

if grad_max < eta_max
    grad_scale=1;
else 
    grad_scale =  grad_max/eta_max;
    
for i_cell = 1:n_cell
    for i_field = 1:length(fldnames)
         new_gradient(i_cell).(fldnames{i_field}).mean= ...
             new_gradient(i_cell).(fldnames{i_field}).mean/grad_scale;
         new_gradient(i_cell).(fldnames{i_field}).sigma= ...
             new_gradient(i_cell).(fldnames{i_field}).sigma/grad_scale;
    end
end
end

