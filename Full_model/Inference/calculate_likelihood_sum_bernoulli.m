function [likelihood] = calculate_likelihood_sum_bernoulli(...
    n_events, probabilities)
%
%      gammas=[gamma0; gamma_sample];
%      probabilities=[1 cells_probabilities(i_trial,:)]';
%t1=toc;
event_prob = probabilities;
stimulated_index = find(event_prob > 1e-4);
n_stimulated =length(stimulated_index);
event_prob = event_prob(stimulated_index);
not_fire_prob = 1-event_prob;

if n_events == 0
    likelihood = prod(not_fire_prob);
elseif n_events > n_stimulated
    % the gamma sample is not feasible
    likelihood= 1e-20;
else
    combinations_of_event_sources = combnk(1:n_stimulated,n_events);
    prob_combs = zeros(size(combinations_of_event_sources,1),1);
    for k = 1:size(combinations_of_event_sources,1)
        cond_prob = not_fire_prob;
        cond_prob(combinations_of_event_sources(k,:))=event_prob(combinations_of_event_sources(k,:));
        prob_combs(k)=prod(cond_prob);
    end
    likelihood = sum(prob_combs);
 end
% %%
% gammas=gamma_sample;probabilities=designs(i_trial,:)';
% n_events=3;
% 
% %%
%  likelihood
% 0.3489+0.4636+0.1689+0.0187
% 


