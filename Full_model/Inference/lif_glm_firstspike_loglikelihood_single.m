function [loss] = lif_glm_firstspike_loglikelihood_single(gain_one,gamma_one,...
    prob_relevant,v_trace_relevant,relevant_trials,mpp_relevant,stim_relevant,...
   v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob)
 
    n_grid = length(v_trace_relevant{1});
   unique_stims = unique(stim_relevant);
     n_unique_stim = length(unique_stims);
     prob_first_spike_delayed_grid=zeros([n_grid,n_unique_stim]);
    for i_stim = 1:n_unique_stim
         i=min(find(stim_relevant==unique_stims(i_stim) ));
         v_trace_one=v_trace_relevant{i}; 
         [prob_first_spike_delayed] = voltage_to_prob(gain_one,  v_trace_one,...
             v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob);
         prob_first_spike_delayed_grid(:,i_stim) =prob_first_spike_delayed;
    end
     
%     figure(1)
% histogram([mpp_relevant.times])
% xlim([0,300])
% figure(2)
% plot(prob_first_spike_delayed_grid(:,1))
    

    loglklh=zeros(length(relevant_trials),1);
    for i =1:length(relevant_trials)
        i_stim = find(unique_stims == stim_relevant(i));
        prob_first_spike_delayed=prob_first_spike_delayed_grid(:,i_stim);
       
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
    
    
    