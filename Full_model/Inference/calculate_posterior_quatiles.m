function [twentyfifth, seventyfifth] = calculate_posterior_quatiles(alpha,beta,lower_bound,upper_bound)
%  alpha=parameter_history.alpha(:,end);
%  beta=parameter_history.beta(:,end);

vf_type =2; %logit normal
if vf_type == 1
    mean_temp= (upper_bound+ (upper_bound-lower_bound)./(1+beta./alpha));
    var_temp=ones(length(alpha),1);
   
elseif vf_type==2
    normal_samples =norminv([0.25; 0.75]);% for evaluating mean & variance gamma 

    % obtain estimates
    twentyfifth=zeros(length(alpha),1);
    seventyfifth=zeros(length(alpha),1);
    
    for i_temp = 1:length(alpha)
        gamma_i=beta(i_temp)*normal_samples+alpha(i_temp);
        gamma_i=exp(gamma_i)./(1+exp(gamma_i));
        gamma_i=(lower_bound+ (upper_bound-lower_bound)*gamma_i);
        
        twentyfifth(i_temp)= gamma_i(1);
        seventyfifth(i_temp)= gamma_i(2);
    end
    
end
end