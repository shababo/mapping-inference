function [prob_relevant] = subtract_one_cell(gain_one_old,gamma_one_old,...
    prob_relevant,v_trace_relevant,relevant_trials,mpp_relevant,...
   v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob,subtract)
     % subtract the cell's effect on the probability
     n_grid = length(v_trace_relevant{1});
    for i = 1:length(relevant_trials)
        prob_firing=zeros([n_grid,1]);
        prob_first_spike=zeros([n_grid,1]);
        prob_first_spike_delayed=zeros([n_grid,1]);
        not_spike_prob=1;
        for i_grid = 2:n_grid
            prob_firing(i_grid)=...
                linkfunc{3}(gain_one_old*v_trace_relevant{i}(i_grid)-...
                v_th_known_one);
            not_spike_prob = not_spike_prob*(1-prob_firing(i_grid -1));
            prob_first_spike(i_grid) =not_spike_prob*prob_firing(i_grid);
        end
        for i_grid = 1:n_grid
            idx_time = max(i_grid-max_delay,1): min(i_grid-min_delay,delay_params.n_grid);
            idx_delay = -( (min(idx_time)-i_grid) : (max(idx_time)-i_grid))+1;
            temp=0;
            for i_time = 1:length(idx_time)
                temp=temp+...
                    prob_first_spike(idx_time(i_time))*delay_prob(idx_delay(i_time));
            end
            prob_first_spike_delayed(i_grid) =temp;
        end
        if subtract
            prob_relevant{i}.no_spike=prob_relevant{i}.no_spike-...
                log(1-gamma_one_old+gamma_one_old*(1-sum(prob_first_spike_delayed)));
            if isempty(mpp_relevant(i).times)==false
                prob_relevant{i}.spikes=prob_relevant{i}.spikes-...
                    gamma_one_old*prob_first_spike_delayed(round(mpp_relevant(i).times));
            end
        else
                prob_relevant{i}.no_spike=prob_relevant{i}.no_spike+...
                log(1-gamma_one_old+gamma_one_old*(1-sum(prob_first_spike_delayed)));
            if isempty(mpp_relevant(i).times)==false
                prob_relevant{i}.spikes=prob_relevant{i}.spikes+...
                    gamma_one_old*prob_first_spike_delayed(round(mpp_relevant(i).times));
            end
       
        end
    end
    
    