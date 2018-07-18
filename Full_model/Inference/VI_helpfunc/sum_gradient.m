function [new_gradient]=sum_gradient(gradients,eta,eta_max,iteration)
%%
step_size =  eta*(iteration)^(-1);
%     step_size =  eta*(iteration)^(-1);
if eta_max>step_size
    eta_max=step_size;
end
%%
fldnames = fieldnames(gradients(1,1));
n_cell = size(gradients,2);
S=size(gradients,1);
new_gradient=struct;
grad_max=zeros(n_cell,1);
for i_cell = 1:n_cell
    for i_field = 1:length(fldnames)
        %         if ~strcmp(fldnames{i_field},'shapes')
        if strcmp(fldnames{i_field},'PR') | strcmp(fldnames{i_field},'gain')
            
            this_field_sigma_f=zeros(S,1);
            this_field_mean_f=zeros(S,1);
            for s = 1:S
                this_field_mean_f(s) = gradients(s,i_cell).(fldnames{i_field}).mean_f;
                %             this_field_mean_h(s) = gradients(s,i_cell).(fldnames{i_field}).mean_h;
                this_field_sigma_f(s) = gradients(s,i_cell).(fldnames{i_field}).sigma_f;
                %             this_field_sigma_h(s) = gradients(s,i_cell).(fldnames{i_field}).sigma_h;
            end
            if isfield(gradients(s,i_cell).(fldnames{i_field}), 'mean_h')
                this_field_sigma_h=zeros(S,1);
                this_field_mean_h=zeros(S,1);
                for s = 1:S
                    this_field_mean_h(s) = gradients(s,i_cell).(fldnames{i_field}).mean_h;
                    this_field_sigma_h(s) = gradients(s,i_cell).(fldnames{i_field}).sigma_h;
                    
                end
            end
            if strcmp(gradients(s,i_cell).(fldnames{i_field}).type, 'spiked-logit-normal')
                this_field_prob_logit_f=zeros(S,1);
                for s = 1:S
                    this_field_prob_logit_f(s) = gradients(s,i_cell).(fldnames{i_field}).prob_logit_f;
                end
            end
            new_gradient(i_cell).(fldnames{i_field})=struct;
            new_gradient(i_cell).(fldnames{i_field}).type=gradients(s,i_cell).(fldnames{i_field}).type;
            if strcmp(gradients(1,i_cell).(fldnames{i_field}).type,'individual')
                if isfield(gradients(s,i_cell).(fldnames{i_field}), 'mean_h')
                    a_mean= quick_cov( this_field_mean_f, this_field_mean_h)/quick_cov( this_field_mean_h, this_field_mean_h);
                    a_sigma=quick_cov( this_field_sigma_f, this_field_sigma_h)/quick_cov( this_field_sigma_h, this_field_sigma_h);
                    mh=mean(this_field_mean_h);sh=mean(this_field_sigma_h);
                else
                    a_mean=0;a_sigma=0;
                    mh=0;sh=0;
                end
                new_gradient(i_cell).(fldnames{i_field}).mean= ...
                    step_size*(mean(this_field_mean_f))-a_mean*mh;
                new_gradient(i_cell).(fldnames{i_field}).sigma= ...
                    step_size*(mean(this_field_sigma_f))-a_sigma*sh;
                if strcmp(gradients(s,i_cell).(fldnames{i_field}).type, 'spiked-logit-normal')
                    new_gradient(i_cell).(fldnames{i_field}).prob_logit= ...
                        step_size*(mean(this_field_prob_logit_f)); %-a*mean(this_field_sigma_h));
                end
            else % if it is a common parameter:
                new_gradient(i_cell).(fldnames{i_field}).mean= ...
                    step_size*(mean(this_field_mean_f));%-a*mean(this_field_mean_h));
                new_gradient(i_cell).(fldnames{i_field}).sigma= ...
                    step_size*(mean(this_field_sigma_f));%-a*mean(this_field_sigma_h));
                if strcmp(gradients(s,i_cell).(fldnames{i_field}).type, 'spiked-logit-normal')
                    new_gradient(i_cell).(fldnames{i_field}).prob_logit= ...
                        step_size*(mean(this_field_prob_logit_f));%-a*mean(this_field_sigma_h));
                end
            end
            grad_max(i_cell)=max([grad_max(i_cell) abs( new_gradient(i_cell).(fldnames{i_field}).mean) abs( new_gradient(i_cell).(fldnames{i_field}).sigma)] );
            if strcmp(gradients(s,i_cell).(fldnames{i_field}).type, 'spiked-logit-normal')
                grad_max(i_cell)=max([grad_max(i_cell) abs( new_gradient(i_cell).(fldnames{i_field}).prob_logit) ]);
            end
        else
            %             n_locs = length(gradients(s,i_cell).(fldnames{i_field}).mean_f);
            %              this_field_sigma_f=zeros(S,n_locs);
            %             this_field_mean_f=zeros(S,n_locs);
            %             for s = 1:S
            %                 for i_loc = 1:n_locs
            %                 this_field_mean_f(s,i_loc) = gradients(s,i_cell).(fldnames{i_field}).mean_f(i_loc);
            %                 %             this_field_mean_h(s) = gradients(s,i_cell).(fldnames{i_field}).mean_h;
            %                 this_field_sigma_f(s,i_loc) = gradients(s,i_cell).(fldnames{i_field}).sigma_f(i_loc);
            %                 %             this_field_sigma_h(s) = gradients(s,i_cell).(fldnames{i_field}).sigma_h;
            %                 end
            %             end
            %             if isfield(gradients(s,i_cell).(fldnames{i_field}), 'mean_h')
            %                 this_field_sigma_h=zeros(S,n_locs);
            %                 this_field_mean_h=zeros(S,n_locs);
            %                 for s = 1:S
            %
            %                 for i_loc = 1:n_locs
            %                     this_field_mean_h(s,i_loc) = gradients(s,i_cell).(fldnames{i_field}).mean_h(i_loc);
            %                     this_field_sigma_h(s,i_loc) = gradients(s,i_cell).(fldnames{i_field}).sigma_h(i_loc);
            %                 end
            %                 end
            %             end
            %             new_gradient(i_cell).(fldnames{i_field})=struct;
            %             new_gradient(i_cell).(fldnames{i_field}).type=gradients(s,i_cell).(fldnames{i_field}).type;
            %             for i_loc = 1:n_locs
            %                 if isfield(gradients(s,i_cell).(fldnames{i_field}), 'mean_h')
            %                     a= quick_cov( this_field_mean_f(:,i_loc), this_field_mean_h(:,i_loc))+quick_cov( this_field_sigma_f(:,i_loc), this_field_sigma_h(:,i_loc));
            %                     a=a/(quick_cov( this_field_mean_h(:,i_loc), this_field_mean_h(:,i_loc))+quick_cov( this_field_sigma_h(:,i_loc), this_field_sigma_h(:,i_loc)));
            %                     mh=mean(this_field_mean_h(:,i_loc));sh=mean(this_field_sigma_h(:,i_loc));
            %                 else
            %                     a=0;
            %                     mh=0;sh=0;
            %                 end
            %                 new_gradient(i_cell).(fldnames{i_field}).mean(i_loc)= ...
            %                     step_size*(mean(this_field_mean_f(:,i_loc)))-a*mh;
            %                 new_gradient(i_cell).(fldnames{i_field}).sigma(i_loc)= ...
            %                     step_size*(mean(this_field_sigma_f(:,i_loc)))-a*sh;
            %             end
            %             grad_max(i_cell)=max([grad_max(i_cell) abs( new_gradient(i_cell).(fldnames{i_field}).mean) abs( new_gradient(i_cell).(fldnames{i_field}).sigma)] );
            
        end
        
        
    end
    if grad_max(i_cell) < eta_max
        grad_scale=1;
    else
        grad_scale =  grad_max(i_cell) /eta_max;
        for i_field = 1:length(fldnames)
            if  strcmp(fldnames{i_field},'PR') | strcmp(fldnames{i_field},'gain')
                new_gradient(i_cell).(fldnames{i_field}).mean= ...
                    new_gradient(i_cell).(fldnames{i_field}).mean/grad_scale;
                new_gradient(i_cell).(fldnames{i_field}).sigma= ...
                    new_gradient(i_cell).(fldnames{i_field}).sigma/grad_scale;
            end
        end
    end
    
end


