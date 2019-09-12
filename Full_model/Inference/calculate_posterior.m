function [post_stat] = calculate_posterior(params,quantile_prob)
%%
% quantile_prob=[0.05 0.5 0.95];
% params=neighbourhood.neurons.params(end);
%%
%% it outputs the median instead of mean
fldnames = fieldnames(params(1));

post_stat=struct;
for i_field = 1:length(fldnames)
    this_params=params.(fldnames{i_field});
    %         this_sample=struct;
    normal_qt =norminv(quantile_prob);%
    normal_samples =normrnd(0,1, [100 1]);% for evaluating mean & variance gamma
% <<<<<<< HEAD
%     if ~strcmp(fldnames{i_field},'shapes') && ~strcmp(fldnames{i_field},'GP_minimal_variance')
% =======
    if ~strcmp(fldnames{i_field},'shapes') & ~strcmp(fldnames{i_field},'z') & ~strcmp(fldnames{i_field},'xy')  
% >>>>>>> 12c2460e7b5cc64bb7eb5fc4ac439dfea46acc2d
        normal_qt=this_params.mean+normal_qt*exp(this_params.log_sigma/2);
        normal_samples=this_params.mean+normal_samples*exp(this_params.log_sigma/2);
        switch this_params.dist
            case 'normal'
                this_qt=normal_qt;
                this_sample=normal_samples;
            case 'log-normal'
                this_qt=exp(normal_qt);
                this_sample=exp(normal_samples);
            case 'logit-normal'
                this_qt=exp(normal_qt)./(1+exp(normal_qt))*...
                    (this_params.bounds.up-this_params.bounds.low)+this_params.bounds.low;
                this_sample=exp(normal_samples)./(1+exp(normal_samples))*...
                    (this_params.bounds.up-this_params.bounds.low)+this_params.bounds.low;
            case 'spiked-logit-normal'
                zero_prob = exp(this_params.prob_logit)./(exp(this_params.prob_logit)+1);
                this_qt=exp(normal_qt)./(1+exp(normal_qt))*...
                    (this_params.bounds.up-this_params.bounds.low)+this_params.bounds.low;
                this_sample=exp(normal_samples)./(1+exp(normal_samples))*...
                    (this_params.bounds.up-this_params.bounds.low)+this_params.bounds.low;
        end
        post_stat.(fldnames{i_field}).mean=this_qt(2);
        post_stat.(fldnames{i_field}).sd=std(this_sample);
        post_stat.(fldnames{i_field}).upper_quantile=this_qt(3);
        post_stat.(fldnames{i_field}).lower_quantile=this_qt(1);
    else %if strcmp(fldnames{i_field},'shapes') % shape parameters:
        post_stat.(fldnames{i_field})=struct;
        if ~isempty(this_params.mean)
            
            if strcmp(this_params.dist,'mvn')
                this_mean=(this_params.bounds.up-this_params.bounds.low).*exp(this_params.mean)./(1+exp(this_params.mean))...
                    +this_params.bounds.low;
                Dmat= diag(exp(-2*this_params.log_sigma) );
                
                Sigma_tilde=inv(this_params.Sigma_tilde+Dmat);
                this_sigma = sqrt(diag(Sigma_tilde));
            end
            for i_loc = 1:length(this_params.mean)
                  normal_qt =norminv(quantile_prob);%
                  normal_samples =normrnd(0,1, [100 1]);% for evaluating mean & variance gamma
                switch this_params.dist
                    case 'mvn'
                        normal_qt=this_mean(i_loc)+normal_qt*this_sigma(i_loc);
                        normal_samples=this_mean(i_loc)+normal_samples*this_sigma(i_loc);
                        this_qt=normal_qt;
                        this_sample=normal_samples;
                     case 'logit-normal'
                         normal_qt=this_params.mean(i_loc)+normal_qt*exp(this_params.log_sigma(i_loc));
                         normal_samples=this_params.mean(i_loc)+normal_samples*exp(this_params.log_sigma(i_loc));
                         this_qt=exp(normal_qt)./(1+exp(normal_qt))*...
                            (this_params.bounds.up(i_loc)-this_params.bounds.low(i_loc))+this_params.bounds.low(i_loc);
                        this_sample=exp(normal_samples)./(1+exp(normal_samples))*...
                            (this_params.bounds.up(i_loc)-this_params.bounds.low(i_loc))+this_params.bounds.low(i_loc);
                end
                post_stat.(fldnames{i_field}).mean(i_loc)=this_qt(2);
                post_stat.(fldnames{i_field}).sd(i_loc)=std(this_sample);
                post_stat.(fldnames{i_field}).upper_quantile(i_loc)=this_qt(3);
                post_stat.(fldnames{i_field}).lower_quantile(i_loc)=this_qt(1);
                post_stat.(fldnames{i_field}).locations(i_loc,:)=this_params.locations(i_loc,:);
            end
        end
    end
end










