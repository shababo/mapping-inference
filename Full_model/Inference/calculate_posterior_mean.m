function [mean_temp, var_temp] = calculate_posterior_mean(alpha,beta,lower_bound,upper_bound)
%  alpha=parameter_history.alpha(:,end);
%  beta=parameter_history.beta(:,end);

vf_type =2; %logit normal
if vf_type == 1
    mean_temp= (upper_bound+ (upper_bound-lower_bound)./(1+beta./alpha));
    var_temp=ones(length(alpha),1);
   
elseif vf_type==2
    normal_samples =normrnd(0,1, [100 1]);% for evaluating mean & variance gamma 

    % obtain estimates
    mean_temp=zeros(length(alpha),1);
    var_temp=zeros(length(alpha),1);
    
    for i_temp = 1:length(alpha)
        gamma_i=beta(i_temp)*normal_samples+alpha(i_temp);
        gamma_i=exp(gamma_i)./(1+exp(gamma_i));
        gamma_i=(lower_bound+ (upper_bound-lower_bound)*gamma_i);
        
        mean_temp(i_temp)= mean(gamma_i);
        var_temp(i_temp)= var(gamma_i);
    end
    
end
end