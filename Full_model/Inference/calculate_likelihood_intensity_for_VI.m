function [loss] = calculate_likelihood_intensity_for_VI(...
    mpp_this_trial, prob_this_trial)
intensity_this_trial=mpp_this_trial.intensity;
intensity_pred = sum(prob_this_trial,1);

loss= sum((intensity_this_trial - intensity_pred).^2);

% if likelihood< 1e-20
%         likelihood=1e-20;
%     end
%  loss=log(likelihood);
end

    
    
    
    