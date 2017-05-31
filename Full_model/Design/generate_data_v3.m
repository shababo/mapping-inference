function [mpp_new, presynaptic_events, background_events] = generate_data_v3(...
    locations_this_batch,powers_this_batch,pi_dense_all,k_minimum,...
    cell_params, shape_template,...
    I_e_vect, data_params,bg_params,trials_specific_variance)

num_this_batch = size(locations_this_batch,1);
% Include power here:
X_next_all = zeros(size(pi_dense_all,1), num_this_batch);
for l = 1:num_this_batch
    for m = 1:size(locations_this_batch,2)
        if isnan(locations_this_batch(l,m))
            
        else
            X_next_all(:,l) = X_next_all(:,l)+( pi_dense_all(:,locations_this_batch(l,m)).*powers_this_batch(l));
        end
       
    end
end
evoked_k = X_next_all;
presynaptic_events = cell(size(evoked_k,2),size(evoked_k,1));
presynaptic_amplitudes = cell(size(evoked_k,2),size(evoked_k,1));
t_vect=0:data_params.dt:data_params.T-1;
    
for i_cell = 1:size(evoked_k,1)

    if cell_params.gamma(i_cell) > 0
        fprintf('%d\n', i_cell);
        mean_amp_this_trial = abs(normrnd(cell_params.amplitudes(i_cell),cell_params.sigma_across(i_cell)));
        sigma=cell_params.sigma_within(i_cell);
        
        params_sim.V_th= cell_params.V_th(i_cell);
        params_sim.V_reset = cell_params.V_reset(i_cell);
        params_sim.g = shape_template( cell_params.shape_gain(i_cell)).g;
        params_sim.gain = shape_template( cell_params.shape_gain(i_cell)).optical_gain;
        funcs.invlink = @invlink_test;%@(resp) log(1 + exp(resp));%@(x) exp(x);
        
        for i_trial = 1:size(evoked_k,2)
            k = evoked_k(i_cell,i_trial);
            stim = I_e_vect*k;
            presynaptic_events{i_trial, i_cell} = [];
            if max(stim)*params_sim.gain > k_minimum
            [V_vect, spikes]  = lif_glm_sim_v2(stim,params_sim,funcs);
            presynaptic_events{i_trial, i_cell} = t_vect(find(spikes));
            presynaptic_amplitudes{i_trial, i_cell} = abs(normrnd(mean_amp_this_trial,sigma,sum(spikes)));
            end
        end
    end
end

% Draw the spontaneous events:
background_events = cell(size(evoked_k,2),1);
mu_bg = 1/bg_params.firing_rate;
for i_trial = 1:size(evoked_k,2)
    background_events{i_trial}=[];
    T_remain =data_params.T;
    R = exprnd(mu_bg);
    while R < T_remain
        background_events{i_trial}=[background_events{i_trial} R];
        T_remain = T_remain - R;
        R = exprnd(mu_bg);
    end
end

%--------------------------------------
% Turn the presynaptic events into postsynaptic events
%   - with synaptic failure
%   - draw event sizes

% sample "ground truth" stimulations
for l = 1:num_this_batch
    mpp_n = struct();
    if isempty(background_events{l}) ==  0
        mpp_n.event_times = background_events{l};
        mpp_n.amplitudes = exp(normrnd(bg_params.mean,bg_params.sigma , ...
            [1 length(background_events{l})]));
        mpp_n.assignments =  zeros(1,length(background_events{l}));
    else
        mpp_n.event_times=[];
        mpp_n.amplitudes=[];
        mpp_n.assignments=[];
    end
    for i_cell = 1:size(evoked_k,1)
        if cell_params.gamma(i_cell)  > 0
            if isempty(presynaptic_events{l,i_cell}) ==  0
                for i = 1:length(presynaptic_events{l,i_cell})
                    
                    if rand(1) < cell_params.gamma(i_cell)
                        mpp_n.event_times = [mpp_n.event_times presynaptic_events{l,i_cell}(i)];
                        mpp_n.amplitudes = [mpp_n.amplitudes presynaptic_amplitudes{l, i_cell}(i)];
                        mpp_n.assignments = [mpp_n.assignments i_cell];
                    end
                end
            end
        end
    end
    if l==1
        mpp_new = mpp_n;
    else
        mpp_new(l) = mpp_n;
    end
end

