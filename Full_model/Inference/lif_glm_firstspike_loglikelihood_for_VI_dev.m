function [loss] = lif_glm_firstspike_loglikelihood_for_VI_dev(...
    mpp_this_trial, prob_this_trial)
            
%             mpp_this_trial=mpp(i_trial);
            
%prob_this_trial
% n_stimulated = prob_this_trial(:,end)>1e-2;
% stimmed_cells=find(sum(prob_this_trial,2)>1e-2);
% prob_this_trial=prob_this_trial(stimmed_cells,:);
n_stimulated= size(prob_this_trial,1);
n_events=length(mpp_this_trial.event_times);
event_prob = prob_this_trial(:,end);
not_fire_prob = 1-event_prob;
%n_grid=size(prob_this_trial,2);
      
if n_events == 0
    likelihood = prod(not_fire_prob);
elseif n_events > n_stimulated
    % the gamma sample is not feasible
    likelihood= 1e-20;
else
    combinations_of_event_sources = combnk(1:n_stimulated,n_events);
    prob_combs = zeros(size(combinations_of_event_sources,1)*factorial(n_events),1);
    
    for k = 1:size(combinations_of_event_sources,1)
        fired_cells=combinations_of_event_sources(k,:);
        orders_of_cells = perms(fired_cells);
        for j= 1:size(orders_of_cells,1)
            cond_prob = not_fire_prob;
            for i_event = 1:n_events
                cond_prob(orders_of_cells(i_event))=prob_this_trial(orders_of_cells(i_event), ...
                   i_event);
            end
            prob_combs( (k-1)*size(orders_of_cells,1)+j)=prod(cond_prob);
        end
    end
    likelihood = sum(prob_combs);
    
end
if likelihood< 1e-20
        likelihood=1e-20;
    end
 loss=log(likelihood);
end

    
    
    
    