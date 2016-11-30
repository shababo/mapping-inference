%% Option II: fit a full model but use the online update 
% This might be slow
% Good news is that it does not depend on Stage I too much 
% To improve the effiiency, we need to reduce the density of the grid 

%% 

% initialize alpha, mu, and s_sq to the estimates from the reduced model
     output= struct([]);
   for j = 1:(num_threshold-1)
      output(j).alpha = .1*ones(K_z+1,1);
        output(j).mu = zeros(K_z+1,1);
        output(j).s_sq = ones(K_z+1,1);
        output(j).w_estimate = zeros(K_z+1,1);
   end   
obj_function = @joint_sparsity_weight_entropy; %
X_g = zeros(N*num_trials_batch,K_z);
Y_g = zeros(N*num_trials_batch,num_threshold-1);
        
tic
locations_record = [];
time_record = zeros(ND,1);
for t = 1:N
    fprintf('Batch %d of %d\n', t, ND);
    % Designing trials
    if t < N_r % The random design
        X_next=[];
        X_next_all=[];
        for l = 1:num_trials_batch
            locations_next=zeros(num_sources,1);
            for l = 1:num_sources
                locations_next(l) = randsample(x_group{l}, 1);
            end
            locations_record = [locations_record locations_next'];
            X_next = [X_next, min(.95,sum(pi_Z_grid(:, locations_next) ,2))];
            X_next_all = [X_next_all, min(.95,sum(pi_Z_all(:,locations_next),2))];
        end
    else % The optiml design
        %--------------------------------------------------%
        % Part I: evaluate the entropies
        delta_H = zeros(K_z,1);
        for j = 1:(num_threshold-1)
            % compute current objective function value
            H_current = per_neuron_joint_entropy(output(j).alpha(2:end), ...
                output(j).mu(2:end), output(j).s_sq(2:end));
            % for each neuron, calculate the expected change in the objective function
            hyperparam_sigma_n = std( sqrt(Y_g(1:(t*num_trials_batch),j)));
            if hyperparam_sigma_n == 0
                hyperparam_sigma_n = 1;
            end
            
            H_expected = approximate_expected_joint_entropy_single_neuron(...
                output(j).alpha(2:end), output(j).mu(2:end), ...
                output(j).s_sq(2:end), hyperparam_sigma_n, num_samples);
            delta_H = delta_H + H_current-H_expected;
            %fprintf('%d %d \n', max(H_current),  max(H_expected));
        end
        entropy_locations = weights_mat_full'*delta_H;
        
        idx = [];
        for i = 1:num_peaks
            [m_id temp_idx] = max(entropy_locations);
            idx = [idx temp_idx];
            entropy_locations = entropy_locations - inner_normalized_products(:,temp_idx)*m_id;
        end
        % Find multiple treaments:
        X_next=[];
        X_next_all=[];
        for l = 1:num_trials_batch
            locations_next = randsample(size(idx,2), num_sources);
            % Add some
            locations_record = [locations_record idx(locations_next)];
            X_next = [X_next, min(.95,sum(pi_dense_grid(:, idx(locations_next)) ,2))];
            X_next_all = [X_next_all, min(.95,sum(pi_dense_all(:,idx(locations_next)),2))];
        end
    end
    
    % update X matrix
    X_g((t-1)*num_trials_batch+(1:num_trials_batch),:)  = X_next';
      
    %--------------------------------------
    % Draw a new sample from 
    pi_k_spike = X_next_all;
    pi_k_spike(pi_k_spike > .65) = 1;
    
    % firing delay means and variances
    d_mean_nk = d_mean0 + (1.5 - X_next_all)*d_mean_coef;
    d_sigma_nk = d_sigma0 + (1 - X_next_all)*d_sigma_coef;
    
    % sample "ground truth" firing delay
    D = exprnd(d_mean_nk/data_params.dt) + evoked_params.stim_start;
    % sample "ground truth" stimulations
    X = rand(size(D,1),size(D,2)) < pi_k_spike; %.2
    X(D > 2000) = 0;
    for l = 1:num_trials_batch
        firing_neurons = X(:,l) & all_amplitudes > 0;
        if any(all_amplitudes(firing_neurons) > 0)
            evoked_params.times = D(firing_neurons);
            evoked_params.a = all_amplitudes(firing_neurons);
            evoked_params.tau_r = all_tau_rise(firing_neurons)/data_params.dt;
            evoked_params.tau_f = all_tau_fall(firing_neurons)/data_params.dt;
        else
            evoked_params.times = [];
            evoked_params.a = [];
            evoked_params.tau_r = [];
            evoked_params.tau_f = [];
        end
        [~, mpp_n] = gen_trace_noise(data_params,bg_params,evoked_params);
        if (t == 1) && (l==1)
            mpp_new = mpp_n;
        else
            mpp_new( (t-1)*num_trials_batch+l  ) = mpp_n;
        end
    end
    
    % Count the events in each amplitude bins:
    related_mpp_n=struct();
    for l = 1:num_trials_batch
        if size(mpp_new( (t-1)*num_trials_batch+l  ).event_times,2) > 0
            indices = mpp_new( (t-1)*num_trials_batch+l  ).event_times>evoked_params.stim_start  & mpp_new( (t-1)*num_trials_batch+l  ).event_times< (400+evoked_params.stim_start);
            related_mpp_n(l).amplitudes = mpp_new( (t-1)*num_trials_batch+l  ).amplitudes(indices);
            related_mpp_n(l).event_times = mpp_new( (t-1)*num_trials_batch+l  ).event_times(indices);
        end
        for j = 1:(num_threshold-1)
            Y_g((t-1)*num_trials_batch+l,j) = sum(related_mpp_n(l).amplitudes>amplitude_threshold(j) & related_mpp_n(l).amplitudes<(amplitude_threshold(j+2)+0.01));
        end
    end
    
    %---------------------------------------------------------%
    % Part III:
    % Fit the VB model to update the parameters:
    % Initialization:
    output_warm = output;
    X_temp = X_g((t-1)*num_trials_batch+(1:num_trials_batch),:);
    X_temp = [ones(size(X_temp,1),1) X_temp];
    for j = 1:size(Y_g,2)
        %fprintf('%d',j);
        Y_n = Y_g((t-1)*num_trials_batch+(1:num_trials_batch),j);
        % Variance stablization transformation 
        Y_n = sqrt(Y_n);
        
        if  sum(Y_n)==0
            hyperparam_sigma_n = 1;
        else
            hyperparam_sigma_n = std(Y_n);
        end
        
        hyperparam_p_connected = [output_warm(j).alpha];
        
        hyperparam_eta =  output_warm(j).mu;
        hyperparam_sigma_s = sqrt(output_warm(j).s_sq);
        
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
        output(j).w_estimate = alpha_tmp.*mu_tmp;
    end
    %-------------------------------------------%
    time_record(t) = toc;
end


