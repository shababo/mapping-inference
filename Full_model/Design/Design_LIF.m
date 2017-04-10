tic;
time_start = toc;
fprintf('Batch %d of %d\n', t, N);
%tic
% Designing trials
if t == 1
    
    num_this_batch = num_trials_first;
    index_this_batch = 1:num_trials_first;
    start_this_batch = 0;
    X_next = zeros(size(pi_dense_local,1), num_this_batch);
    X_next_all = zeros(size(pi_dense_all,1), num_this_batch);
    counts_freq =  zeros(size(entropy_locations,1),1);
    
    for l = 1:num_this_batch
        locations_next = randsample(grid_index, num_sources);
        locations_trials(start_this_batch+l,:)=locations_next;
        X_next(:,l) =  sum(pi_dense_local(:, locations_next) ,2);
        X_next_all(:,l) = sum(pi_dense_all(:,locations_next),2);
        counts_freq(locations_next) = counts_freq(locations_next)+1;
    end
        
else
    num_this_batch = num_trials_batch;
    index_this_batch = num_trials_first + (t-2)*num_trials_batch + (1:num_trials_batch);
    start_this_batch = num_trials_first + (t-2)*num_trials_batch;
    
    if design == 0
        % Random design
        %entropy_locations = pi_dense_local'*ones(n_cell_local,1);
        idx = grid_index;
    elseif design == 1
        % Optiml design
        %--------------------------------------------------%
        % Part I: evaluate the entropies
        delta_H = zeros(n_cell_local+1,1);
        for j = 1:(num_threshold)
            H_current = per_neuron_joint_entropy(output(j).alpha, ...
                output(j).mu, output(j).s_sq);
            Y_n = Y_g(1:start_this_batch,j);
            if sqrt_transform
                % Variance stablization transformation
                Y_n = sqrt(Y_n);
            end
            if  sum(Y_n)==0 % prevent zero standard deviation
                hyperparam_sigma_n = 1;
            else
                hyperparam_sigma_n = std(Y_n);
            end
            
            H_expected = approximate_expected_joint_entropy_single_neuron_poisson(...
                output(j).alpha, output(j).mu, ...
                output(j).s_sq, hyperparam_sigma_n, num_samples);
            delta_H = delta_H + H_current-H_expected;   
        end
        entropy_locations = pi_dense_local'*delta_H(2:end);
        
        % Encourage stimulating rare locations (e.g. nearby cells)
        weights_freq= 1.1-counts_freq/max(counts_freq); 
        
        % use 1.1 so that the penalty will not be 0. 
        entropy_locations = weights_freq.*entropy_locations;
        idx = [];
        for i = 1:num_peaks
            [m_id temp_idx] = max(entropy_locations);
            idx = [idx temp_idx];
            entropy_locations = entropy_locations - inner_normalized_products(:,temp_idx)*m_id;
        end
    end
    
    % Find multiple treaments:
    % Permute the index vector
    num_rep = ceil(num_sources*num_this_batch/num_peaks);
    idx_seq = zeros(num_rep*num_peaks,1);
    for l = 1:num_rep
        idx_seq((l-1)*num_peaks+(1:num_peaks)) = randsample(idx, num_peaks);
    end
    X_next = zeros(size(pi_dense_local,1), num_this_batch);
    X_next_all = zeros(size(pi_dense_all,1), num_this_batch);
    for l = 1:num_this_batch
        locations_next = idx_seq((l-1)*num_sources + (1:num_sources));
        locations_trials(start_this_batch+l,:)=locations_next;
        X_next(:,l) = sum(pi_dense_local(:, locations_next) ,2);
        X_next_all(:,l) = sum(pi_dense_all(:,locations_next),2);
        counts_freq(locations_next) = counts_freq(locations_next)+1;
    end
end
    
%--------------------------------------
% Draw new samples from LIF
% Spatial marks
evoked_k = k_offset.* X_next_all;
presynaptic_events = cell(size(evoked_k,2),size(evoked_k,1));
presynaptic_amplitudes = cell(size(evoked_k,2),size(evoked_k,1));
for i_cell = 1:size(evoked_k,1)
    if all_amplitudes(i_cell) > 0
        V_th = all_V_th(i_cell);
        V_reset = all_V_reset(i_cell);
        sigma=all_sigma(i_cell);
        for i_trial = 1:size(evoked_k,2)
            k = evoked_k(i_cell,i_trial);
            presynaptic_events{i_trial, i_cell} = [];
            if k > k_minimum
                %%%DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS
                %spTrain=zeros(t_end,length(I_Stim_vect));% The output spike train
                i=1; %index denoting which element of V is being assigned
                V_vect(i)=0; %first element of V, i.e. value of V at t=0
                %%%%chi-sq shape current
                I_e_vect=[0;I_e(:,num_I_Stim)];
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
                            abs(normrnd(all_amplitudes(i_cell),sigma))];
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
                    if rand(1) < all_gamma(i_cell)
                        mpp_n.event_times = [mpp_n.event_times presynaptic_events{l,i_cell}(i)];
                        mpp_n.amplitudes = [mpp_n.amplitudes presynaptic_amplitudes{l, i_cell}(i)];
                        mpp_n.assignments = [mpp_n.assignments i_cell];
                    end
                end
            end
        end
    end
    
    if (t == 1) && (l==1)
        mpp_new = mpp_n;
    else
        mpp_new( start_this_batch+l  ) = mpp_n;
    end
end

%total_events = [mpp_new.assignments];
%precells = local_index(local_connected);

% Count the events in each amplitude bins:
related_mpp_n=struct();
for l = 1:num_this_batch
    if size(mpp_new( start_this_batch+l  ).event_times,2) > 0
        indices = mpp_new(start_this_batch+l).event_times>evoked_params.stim_start & ...
            mpp_new(start_this_batch+l).event_times< evoked_params.stim_end;
        related_mpp_n(l).amplitudes = mpp_new( start_this_batch+l  ).amplitudes(indices);
        related_mpp_n(l).event_times = mpp_new(start_this_batch+l  ).event_times(indices)-evoked_params.stim_start;
    else
        related_mpp_n(l).amplitudes = [];
        related_mpp_n(l).event_times = [];
    end
end
if t == 1
    quans =  (1/num_threshold)*[0:num_threshold];
    quans(end) = quans(end)-0.01; % prevent heavy-tail
    if mark == 0
        threshold=  quantile([related_mpp_n.amplitudes],quans);
        % Alternatively, use an evenly-spaced grid 
        %threshold= linspace(min([related_mpp_n.amplitudes]) ,max([related_mpp_n.amplitudes]),num_threshold+1);
    elseif mark == 1
        threshold=  quantile([related_mpp_n.event_times], quans);
    end
end
for l = 1:num_this_batch
    for j = 1:num_threshold
        if mark == 0
            Y_g(start_this_batch+l,j) = sum(related_mpp_n(l).amplitudes>threshold(j) & related_mpp_n(l).amplitudes<(threshold(j+1)+0.01));
        elseif mark==1
            Y_g(start_this_batch+l,j) = sum(related_mpp_n(l).event_times >threshold(j) & related_mpp_n(l).event_times <(threshold(j+1)+0.01));
        end
    end
end

X_g(index_this_batch,:) =X_next';

%---------------------------------------------------------%
% Part III:
% Fit the VB model to update the parameters:
% Initialization:
output_warm = output;
 X_temp = k_offset.*X_g(max(1,max(index_this_batch)-length_memory):max(index_this_batch),:);
   X_temp = [ones(size(X_temp,1),1) X_temp];
%X_temp = k_offset.*X_g( 1:max(index_this_batch),:);
%   X_temp = [ones(size(X_temp,1),1) X_temp];

for j = 1:size(Y_g,2)
    %fprintf('%d',j);
   Y_n = Y_g(max(1,max(index_this_batch)-length_memory):max(index_this_batch),j);
    
    %Y_n = Y_g(1:max(index_this_batch),j);
    if sqrt_transform
        % Variance stablization transformation
        Y_n = sqrt(Y_n);
    end
    if  sum(Y_n)==0 % prevent zero standard deviation
        hyperparam_sigma_n = 1;
    else
        hyperparam_sigma_n = std(Y_n);
    end
    
    hyperparam_p_connected = output_ini(j).alpha; % No initial values for alpha
    hyperparam_eta =  output_ini(j).mu;
    hyperparam_sigma_s = ones(n_cell_local+1,1);
    %sqrt(output_warm(j).s_sq);
    %
    options = struct();
    options.alpha =  output_warm(j).alpha;
    options.mu= output_warm(j).mu;
    options.verbose= false;
    options.center= 0;
    
    [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs_general(X_temp, Y_n, ...
        hyperparam_sigma_n, hyperparam_sigma_s, hyperparam_p_connected, ...
        hyperparam_eta, options);
    output(j).alpha = alpha_tmp;
    output(j).mu = mu_tmp;
    output(j).s_sq = s_sq_tmp;
end
%toc
%-------------------------------------------%
time_end = toc;
time_batch = time_end - time_start;
time_record(t) = time_batch;

