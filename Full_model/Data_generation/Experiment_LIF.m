%% Run Batches t_start to t_end
% Need to turn this into a function
% Inputs: output, t_start, t_end, pi_dense_all, num_trials_first, num_trials_batch
% pi_dense_local, num_sources, num_this_batch, num_peaks, threshold
% evoked_params, d_mean_nk, d_mean0, d_mean_coef, d_sigma_nk, d_sigma0, d_sigma_coef

% Generate data using the LIF model 


presynaptic_events = cell(n_trial, n_cell);
presynaptic_amplitudes = cell(n_trial, n_cell);
voltage_traces = cell(n_trial, n_cell);
    
% One cell:
for i_cell = 1:n_cell
    
    
    V_th = all_V_th(i_cell);
    V_reset = all_V_reset(i_cell);
    E_L=all_E_L(i_cell); %resting membrane potential [mV]
    num_I_Stim=1;
    for i_trial = 1:n_trial
        presynaptic_events{i_trial, i_cell} = [];
        presynaptic_amplitudes{i_trial,i_cell} = [];
        k = stimuli_size(i_trial, i_cell);
        if k > 0.01
            %%%DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS
            t_vect=0:data_params.dt:data_params.T;
            V_vect=zeros(1,length(t_vect));
            
            %spTrain=zeros(t_end,length(I_Stim_vect));% The output spike train
            i=1; %index denoting which element of V is being assigned
            V_vect(i)=E_L; %first element of V, i.e. value of V at t=0
            %%%%chi-sq shape current
            I_e_vect=[0;I_e(:,num_I_Stim)];
            for ti=data_params.dt:data_params.dt:data_params.T %loop through values of t in steps of df ms
                V_vect(i+1) = V_vect(i) + ((E_L-V_vect(i))*g + I_e_vect(i)*k)*data_params.dt + ...
                    sqrt(data_params.dt)*normrnd(stoc_mu,stoc_sigma);
                
                if exact_crossing == 1
                    % Modify the following statement to make sure that
                    %if statement below says what to do if voltage crosses threshold
                    if (V_vect(i+1)>V_th) %cell spiked
                        V_vect(i+1)=V_reset; %set voltage back to V_reset
                        presynaptic_events{i_trial, i_cell} = [presynaptic_events{i_trial, i_cell} t_vect(i)];
                        presynaptic_amplitudes{i_trial, i_cell} = [presynaptic_amplitudes{i_trial, i_cell} normrnd(all_amplitudes(i_cell),evoked_params.sigma_a)];
                        %t
                    end  
                else
                    prob_spike = min(1,max(0,V_vect(i+1) - V_th)*data_params.dt);
                    spike_indicator = binornd(1,prob_spike);
                    if spike_indicator %cell spiked
                        V_vect(i+1)=V_reset; %set voltage back to V_reset
                        presynaptic_events{i_trial, i_cell} = [presynaptic_events{i_trial, i_cell} t_vect(i)];
                        presynaptic_amplitudes{i_trial, i_cell} = [presynaptic_amplitudes{i_trial, i_cell} normrnd(all_amplitudes(i_cell),evoked_params.sigma_a)];
                        %t
                    end
                end
                
                i=i+1; %add 1 to index,corresponding to moving forward 1 time step
            end
            voltage_traces{i_trial,i_cell} = V_vect;
        end
        
    end
    %presynaptic_events{i_trial, i_cell};
    fprintf('%d\n', i_cell);
end

%%
% Draw the spontaneous events:

bg_params.firing_rate = 10; % was 50

background_events = cell(n_trial,1);
mu_bg = 1000/bg_params.firing_rate;
for i_trial = 1:n_trial
    background_events{i_trial}=[];
    T_remain =data_params.T;
    R = exprnd(mu_bg);
    while R < T_remain
        background_events{i_trial}=[background_events{i_trial} R];
        T_remain = T_remain - R;
        R = exprnd(mu_bg);
    end
end
%background_events{i_trial};
fprintf('bg\n');
%%
%--------------------------------------
% Turn the presynaptic events into postsynaptic events
%   - with synaptic failure
%   - draw event sizes

for l = 1:n_trial
    mpp_n = struct();
    if isempty(background_events{l}) ==  0
        mpp_n.event_times = background_events{l};
        mpp_n.amplitudes = rand([1 length(background_events{l})])...
            .*(bg_params.a_max -bg_params.a_min)+bg_params.a_min;
        mpp_n.assignments =  zeros(1,length(background_events{l}));
        
    else
        mpp_n.event_times=[];
        mpp_n.amplitudes=[];
        mpp_n.assignments=[];
        
    end
    for i_cell = 1:n_cell
        if all_connected(i_cell) > 0
            if isempty(presynaptic_events{l,i_cell}) ==  0
                for i = 1:length(presynaptic_events{l,i_cell})
                    if rand(1) > evoked_params.failure_prob
                        mpp_n.event_times = [mpp_n.event_times presynaptic_events{l,i_cell}(i)];
                        mpp_n.amplitudes = [mpp_n.amplitudes presynaptic_amplitudes{l, i_cell}(i)];
                        mpp_n.assignments = [mpp_n.assignments i_cell];
                    end
                end
            end
        end
    end
    
    if (l==1)
        mpp_new = mpp_n;
    else
        mpp_new(l) = mpp_n;
    end
    fprintf('%d\n', l);
end
