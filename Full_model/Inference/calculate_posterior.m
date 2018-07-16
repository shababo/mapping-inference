function [post_stat] = calculate_posterior(params,quantile_prob)

%% it outputs the median instead of mean!
fldnames = fieldnames(params(1));

post_stat=struct;
for i_field = 1:length(fldnames)
    this_params=params.(fldnames{i_field});
    %         this_sample=struct;
    normal_qt =norminv(quantile_prob);%
    normal_samples =normrnd(0,1, [100 1]);% for evaluating mean & variance gamma
    if ~strcmp(fldnames{i_field},'shapes')
        normal_qt=this_params.mean+normal_qt*exp(this_params.log_sigma);
        normal_samples=this_params.mean+normal_samples*exp(this_params.log_sigma);
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
        post_stat.(fldnames{i_field}).variance=var(this_sample);
        post_stat.(fldnames{i_field}).upper_quantile=this_qt(3);
        post_stat.(fldnames{i_field}).lower_quantile=this_qt(1);
    else
        post_stat.(fldnames{i_field})=struct;
        if ~isempty(this_params.mean)
            for i_loc = 1:length(this_params.mean)
                normal_qt=this_params.mean(i_loc)+normal_qt*exp(this_params.log_sigma(i_loc));
                normal_samples=this_params.mean(i_loc)+normal_samples*exp(this_params.log_sigma(i_loc));
                switch this_params.dist
                    case 'normal'
                        this_qt=normal_qt;
                        this_sample=normal_samples;
                     case 'logit-normal'
                        this_qt=exp(normal_qt)./(1+exp(normal_qt))*...
                            (this_params.bounds.up(i_loc)-this_params.bounds.low(i_loc))+this_params.bounds.low(i_loc);
                        this_sample=exp(normal_samples)./(1+exp(normal_samples))*...
                            (this_params.bounds.up(i_loc)-this_params.bounds.low(i_loc))+this_params.bounds.low(i_loc);
                end
                post_stat.(fldnames{i_field}).mean(i_loc)=this_qt(2);
                post_stat.(fldnames{i_field}).variance(i_loc)=var(this_sample);
                post_stat.(fldnames{i_field}).upper_quantile(i_loc)=this_qt(3);
                post_stat.(fldnames{i_field}).lower_quantile(i_loc)=this_qt(1);
                post_stat.(fldnames{i_field}).locations(i_loc,:)=this_params.locations(i_loc,:);
            end
        end
    end
end










