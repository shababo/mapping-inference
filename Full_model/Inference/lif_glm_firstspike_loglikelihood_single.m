function [loss] = lif_glm_firstspike_loglikelihood_single(gamma_one,...
     prob_trace_one,stim_unique,stim_relevant,n_round_digit,...
    prob_relevant,relevant_trials,mpp_relevant,...
   v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob)
 
  
%     figure(1)
% histogram([mpp_relevant.times])
% xlim([0,300])
% figure(2)
% plot(prob_first_spike_delayed_grid(:,1))
    

    loglklh=zeros(length(relevant_trials),1);
    for i =1:length(relevant_trials)
        i_stim = find(stim_unique == round(stim_relevant(i),n_round_digit));
        prob_first_spike_delayed=prob_trace_one{i_stim};
       
        if isempty(mpp_relevant(i).times)==true
            temp_prob=...
                1-gamma_one+gamma_one*(1-sum(prob_first_spike_delayed));
            if temp_prob < 1e-8
                temp_prob=1e-8; % avoid singularity
            end
            loglklh(i)= prob_relevant{i}.no_spike+log(temp_prob);
        end %otherwise it's zero
        if isempty(mpp_relevant(i).times)==false
            temp = prob_relevant{i}.spikes+ ...
                gamma_one*prob_first_spike_delayed(round(mpp_relevant(i).times));
            loglklh(i)= sum(log(temp));
        end
    end
    
    loss= sum(loglklh);
    
    
    