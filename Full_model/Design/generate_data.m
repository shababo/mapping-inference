%--------------------------------------
% Draw new samples from LIF
% Spatial marks
function [mpp_new, presynaptic_events, background_events] = generate_data(...
    locations_this_batch,pi_dense_all,k_offset,k_minimum,all_V_th,all_V_reset,all_gamma,...
    all_amplitudes,all_sigma_across,all_sigma_within, ...
    I_e_vect,g,stoc_mu,stoc_sigma, data_params,bg_params,trials_specific_variance)

num_this_batch = size(locations_this_batch,1);
X_next_all = zeros(size(pi_dense_all,1), num_this_batch);
for l = 1:num_this_batch
    X_next_all(:,l) = sum(pi_dense_all(:,locations_this_batch(l,:)),2);
end
evoked_k = k_offset.* X_next_all;
presynaptic_events = cell(size(evoked_k,2),size(evoked_k,1));
presynaptic_amplitudes = cell(size(evoked_k,2),size(evoked_k,1));
t_vect=0:data_params.dt:data_params.T;
    
for i_cell = 1:size(evoked_k,1)
    if all_amplitudes(i_cell) > 0
        V_th = all_V_th(i_cell);
        V_reset = all_V_reset(i_cell);
        mean_amp_this_trial = abs(normrnd(all_amplitudes(i_cell),all_sigma_across(i_cell)));
        sigma=all_sigma_within(i_cell);
        for i_trial = 1:size(evoked_k,2)
            V_vect=zeros(1,length(t_vect));
    
            k = evoked_k(i_cell,i_trial);
            presynaptic_events{i_trial, i_cell} = [];
            if k > k_minimum
                %%%DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS
                %spTrain=zeros(t_end,length(I_Stim_vect));% The output spike train
                i=1; %index denoting which element of V is being assigned
                V_vect(i)=0; %first element of V, i.e. value of V at t=0
                %%%%chi-sq shape current
                for ti=data_params.dt:data_params.dt:data_params.T %loop through values of t in steps of df ms
                    V_vect(i+1) = V_vect(i) + ((0-V_vect(i))*g + I_e_vect(i)*k)*data_params.dt + ...
                        sqrt(data_params.dt)*normrnd(stoc_mu,stoc_sigma);
                    %if statement below says what to do if voltage crosses threshold
                    prob_spike = min(1,max(0,V_vect(i+1) - V_th)*data_params.dt);
                    spike_indicator = binornd(1,prob_spike);
                    if spike_indicator %cell spiked
                        V_vect(i+1)=V_reset; %set voltage back to V_reset
                        presynaptic_events{i_trial, i_cell} = [presynaptic_events{i_trial, i_cell} t_vect(i)];
                        presynaptic_amplitudes{i_trial, i_cell} = [presynaptic_amplitudes{i_trial, i_cell} ...
                            abs(normrnd(mean_amp_this_trial,sigma))];
                    end
                    i=i+1; %add 1 to index,corresponding to moving forward 1 time step
                end
            end
        end
    end
end

% Draw the spontaneous events:
background_events = cell(size(evoked_k,2),1);
mu_bg = 1000/bg_params.firing_rate;
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
        if all_amplitudes(i_cell) > 0
            if isempty(presynaptic_events{l,i_cell}) ==  0
                for i = 1:length(presynaptic_events{l,i_cell})
                    if trials_specific_variance == 1
                        gamma_this_event = exp(presynaptic_amplitudes{l, i_cell}(i)/2)./...
                            (exp(presynaptic_amplitudes{l, i_cell}(i)/2)+1);
                    else
                        gamma_this_event = all_gamma(i_cell);
                    end
                    if rand(1) < gamma_this_event
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

