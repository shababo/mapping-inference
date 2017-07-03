function [loss] = lif_glm_firstspike_loglikelihood_single(gain_one,gamma_one,...
    prob_relevant,v_trace_relevant,relevant_trials,mpp_relevant,...
   v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob)
     n_grid = length(v_trace_relevant{1});

  loglklh=zeros(length(relevant_trials),1);
    for i =1:length(relevant_trials)
        prob_firing=zeros([n_grid,1]);
        prob_first_spike=zeros([n_grid,1]);
        prob_first_spike_delayed=zeros([n_grid,1]);
        not_spike_prob=1;
        for i_grid = 2:n_grid
            prob_firing(i_grid)=...
                linkfunc{3}(gain_one*v_trace_relevant{i}(i_grid)-...
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
        
        if isempty(mpp_relevant(i).times)==true
            temp_prob=...
                1-gamma_one+gamma_one*(1-sum(prob_first_spike_delayed));
            if temp_prob <1e-8
                temp_prob=1e-8; % avoid singularity
            end
            loglklh(i)= prob_relevant{i}.no_spike+log(temp_prob);
        end %otherwise it's zero
        if isempty(mpp_relevant(i).times)==false
            temp = prob_relevant{i}.spikes+ ...
                gamma_one*prob_first_spike_delayed(round(mpp_relevant(i).times));
            loglklh(i)= ...
                loglklh(i)+sum(log(temp));
        end
    end
    loss= exp(sum(loglklh));
    
    