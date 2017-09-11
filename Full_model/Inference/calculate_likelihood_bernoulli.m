function [likelihood] = calculate_likelihood_bernoulli(...
    n_events, gammas, probabilities)

event_prob = gammas.*probabilities;
stimulated_index = find(event_prob > 1e-4);
n_stimulated =length(stimulated_index);
event_prob = event_prob(stimulated_index);
not_fire_prob = 1-event_prob;
if n_events == 0
    likelihood = prod(not_fire_prob); 
else
    likelihood = 1-prod(not_fire_prob); 
end
end

    
    
    
    