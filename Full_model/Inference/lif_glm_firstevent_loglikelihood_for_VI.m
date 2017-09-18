function [loss] = lif_glm_firstevent_loglikelihood_for_VI(...
    mpp_this_trial,...
    prob_this_trial)

%             mpp_this_trial=mpp(i_trial);

%prob_this_trial
n_stimulated = sum( sum(prob_this_trial,2)>1e-2);
stimmed_cells=find(sum(prob_this_trial,2)>1e-2);

prob_this_trial=prob_this_trial(stimmed_cells,:);
n_grid=size(prob_this_trial,2);

if isempty(mpp_this_trial.times)
    event_prob = sum(prob_this_trial,2);
    not_fire_prob = 1-event_prob;
    
    likelihood = prod(not_fire_prob);
else
    not_fire_prob = 1-sum(prob_this_trial(:,1:round(mpp_this_trial.times(1))),2);
    prob_combs=zeros(n_stimulated,1);
    for j = 1:n_stimulated
        cond_prob = not_fire_prob;
        cond_prob(j)=prob_this_trial(j,round(mpp_this_trial.times(1)));
        prob_combs(j)=prod(cond_prob);
    end
    likelihood = sum(prob_combs);
    if likelihood< 1e-20
        likelihood=1e-20;
    end
end
loss=log(likelihood);
end




