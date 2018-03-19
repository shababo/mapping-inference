function [loss] = lif_glm_firstevent_loglikelihood_for_VI(...
    mpp_this_trial,...
    prob_this_trial)
n_stimulated= size(prob_this_trial,1);
n_events=length(mpp_this_trial.event_times);
event_prob = prob_this_trial(:,end);
not_fire_prob = 1-event_prob;


if n_events == 0
    likelihood = prod(not_fire_prob);
else
    prob_combs=zeros(n_stimulated,1);
    for j = 1:n_stimulated
        cond_prob = not_fire_prob;
        cond_prob(j)=prob_this_trial(j,1);
        prob_combs(j)=prod(cond_prob);
    end
    likelihood = sum(prob_combs);
    
end

if likelihood< 1e-20
        likelihood=1e-20;
    end
 loss=log(likelihood);
end

    

