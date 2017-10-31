function [loss] = calculate_loglikelihood_bernoulli(...
    mpp_this_trial, probabilities)
n_events=length(mpp_this_trial.times);
event_prob = sum(probabilities,2);
stimulated_index = find(event_prob > 1e-4);
n_stimulated =length(stimulated_index);
event_prob = event_prob(stimulated_index);
not_fire_prob = 1-event_prob;
if n_events == 0
    likelihood = prod(not_fire_prob); 
else
    likelihood = 1-prod(not_fire_prob); 
end
if likelihood < 1e-6
    likelihood=1e-6;
end
loss=log(likelihood);
end

