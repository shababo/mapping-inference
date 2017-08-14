function [prob_relevant] = subtract_one_cell(gamma_one_old,...
     prob_trace_one,stim_unique,stim_relevant,n_round_digit,...
    prob_relevant,relevant_trials,mpp_relevant,...
   v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob,subtract)

     % subtract the cell's effect on the probability
     % use the fact that the stimuli sizes are categorical 
%      n_grid = length(v_trace_relevant{1});
%      unique_stims = unique(stim_relevant);
%      n_unique_stim = length(unique_stims);
%      prob_first_spike_delayed_grid=zeros([n_grid,n_unique_stim]);
%      
%      for i_stim = 1:n_unique_stim
%          i=min(find(stim_relevant==unique_stims(i_stim) ));
%          v_trace_one=v_trace_relevant{i}; gain_one=gain_one_old;
%          [prob_first_spike_delayed] = voltage_to_prob(gain_one,  v_trace_one,...
%              v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob);
%          prob_first_spike_delayed_grid(:,i_stim) =prob_first_spike_delayed;
%      end
     
    for i = 1:length(relevant_trials)
%         i_stim = find(unique_stims == stim_relevant(i));
%         prob_first_spike_delayed=prob_first_spike_delayed_grid(:,i_stim);
        i_stim = find(stim_unique == round(stim_relevant(i),n_round_digit));
        prob_first_spike_delayed=prob_trace_one{i_stim};
        if subtract
            temp_prob=...
                1-gamma_one_old+gamma_one_old*(1-sum(prob_first_spike_delayed));
            if temp_prob < 1e-8
                temp_prob=1e-8; % avoid singularity
            end
            prob_relevant{i}.no_spike=prob_relevant{i}.no_spike-...
                log(temp_prob);
            if isempty(mpp_relevant(i).times)==false
                prob_relevant{i}.spikes=prob_relevant{i}.spikes-...
                    gamma_one_old*prob_first_spike_delayed(round(mpp_relevant(i).times));
            end
        else
            temp_prob=...
                1-gamma_one_old+gamma_one_old*(1-sum(prob_first_spike_delayed));
            if temp_prob < 1e-8
                temp_prob=1e-8; % avoid singularity
            end
                prob_relevant{i}.no_spike=prob_relevant{i}.no_spike+...
                log(temp_prob);
            if isempty(mpp_relevant(i).times)==false
                prob_relevant{i}.spikes=prob_relevant{i}.spikes+...
                    gamma_one_old*prob_first_spike_delayed(round(mpp_relevant(i).times));
            end
        end
    end
    
    