%% Option II: fit a full model but use the online update
% This might be slow
% Good news is that it does not depend on Stage I too much
% To improve the effiiency, we need to reduce the density of the grid
%%
% initialize alpha, mu, and s_sq to the estimates from the reduced model
for t = t_start:t_end
    time_record(2*t) = toc;
    fprintf('Batch %d of %d\n', t, N);
    %tic
    % Designing trials
    if t == 1
        num_this_batch = num_trials_first;
        index_this_batch = 1:num_trials_first;
        start_this_batch = 0;
        idx_seq = 1:size(Z,1);
        X_next = zeros(size(pi_dense_local,1), num_this_batch);
        X_next_all = zeros(size(pi_dense_all,1), num_this_batch);
        for l = 1:num_this_batch
            locations_next = randsample(idx_seq, num_sources);
            locations_trials(start_this_batch+l,:) = locations_next;
            X_next(:,l) = min(.95,sum(pi_dense_local(:, locations_next) ,2));
            X_next_all(:,l) = min(.95,sum(pi_dense_all(:,locations_next),2));
        end
        X_g(index_this_batch,:)  = X_next';
    else
        num_this_batch = num_trials_batch;
        index_this_batch = num_trials_first + (t-2)*num_trials_batch + (1:num_trials_batch);
        start_this_batch = num_trials_first + (t-2)*num_trials_batch;
        
        if design == 0
            % Random design
            entropy_locations = weights_mat_full'*ones(K_z,1);
            idx = randsample(grid_index, num_peaks);
            %             for i = 1:num_peaks
            %                 temp_idx = randsample(grid_index(entropy_locations>0.5) ,1);
            %                 m_id = entropy_locations(temp_idx);
            %                 idx = [idx temp_idx];
            %                 entropy_locations = entropy_locations - inner_normalized_products(:,temp_idx)*m_id;
            %             end
        elseif design == 1
            % The optiml design
            %--------------------------------------------------%
            % Part I: evaluate the entropies
            delta_H = zeros(K_z+1,1);
            for j = 1:(num_threshold)
                % compute current objective function value
                H_current = per_neuron_joint_entropy(output(j).alpha, ...
                    output(j).mu, output(j).s_sq);
                % for each neuron, calculate the expected change in the objective function
                hyperparam_sigma_n = std( sqrt(Y_g(1:(t*num_trials_batch),j)));
                if hyperparam_sigma_n == 0
                    hyperparam_sigma_n = 1;
                end
                H_expected = approximate_expected_joint_entropy_single_neuron_poisson(...
                    output(j).alpha, output(j).mu, ...
                    output(j).s_sq, hyperparam_sigma_n, num_samples);
                delta_H = delta_H + H_current-H_expected;
                %fprintf('%d %d \n', max(H_current),  max(H_expected));
            end
            entropy_locations = weights_mat_full'*delta_H(2:end);
            
            idx = [];
            for i = 1:num_peaks
                [m_id temp_idx] = max(entropy_locations);
                idx = [idx temp_idx];
                entropy_locations = entropy_locations - inner_normalized_products(:,temp_idx)*m_id;
            end
        else
            % The othogonal design
            [~, ~, ranks] = unique(diagsq(X_g));
            delta_H= 1-2*sigmf(ranks,[sigmoid_a mean(ranks)]);
             
            entropy_locations = weights_mat_full'*delta_H;
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
            % Add some
            %locations_record = [locations_record locations_next];
            locations_trials(start_this_batch+l,:) = locations_next;
            X_next(:,l) = min(.95,sum(pi_dense_local(:, locations_next) ,2));
            X_next_all(:,l) = min(.95,sum(pi_dense_all(:,locations_next),2));
        end
        X_g(index_this_batch,:)  = X_next';
    end
    
    %toc
    
    %--------------------------------------
    % Draw a new sample from
    pi_k_spike = X_next_all;
    
    % firing delay means and variances
    d_mean_nk = d_mean0 + (1.5 - X_next_all)*d_mean_coef;
    d_sigma_nk = d_sigma0 + (1 - X_next_all)*d_sigma_coef;
    
    % sample "ground truth" firing delay
    D = exprnd(d_mean_nk/data_params.dt) + evoked_params.stim_start;
    % sample "ground truth" stimulations
    X = rand(size(D,1),size(D,2)) < pi_k_spike; %.2
    X(D > 2000) = 0;
    for l = 1:num_this_batch
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
            mpp_new( start_this_batch+l  ) = mpp_n;
        end
    end
    %histogram([mpp_new.amplitudes])
    if t == 1
        amplitude_threshold=  quantile([mpp_new.amplitudes], (1/num_threshold)*[0:num_threshold]);
    end
    % Count the events in each amplitude bins:
    related_mpp_n=struct();
    for l = 1:num_this_batch
        if size(mpp_new( start_this_batch+l  ).event_times,2) > 0
            indices = mpp_new(start_this_batch+l).event_times>evoked_params.stim_start ...
                & mpp_new( start_this_batch+l  ).event_times< (400+evoked_params.stim_start);
            related_mpp_n(l).amplitudes = mpp_new( start_this_batch+l  ).amplitudes(indices);
            related_mpp_n(l).event_times = mpp_new(start_this_batch+l  ).event_times(indices);
        else
            related_mpp_n(l).amplitudes = [];
            related_mpp_n(l).event_times = [];
        end
        for j = 1:num_threshold
            Y_g(start_this_batch+l,j) = sum(related_mpp_n(l).amplitudes>amplitude_threshold(j) & related_mpp_n(l).amplitudes<(amplitude_threshold(j+1)+0.01));
        end
    end
    %toc
    %---------------------------------------------------------%
    % Part III:
    % Fit the VB model to update the parameters:
    % Initialization:
    output_warm = output;
    X_temp = X_g(index_this_batch,:);
    X_temp = [ones(size(X_temp,1),1) X_temp];
    for j = 1:size(Y_g,2)
        %fprintf('%d',j);
        Y_n = Y_g(index_this_batch,j);
        if sqrt_transform
            % Variance stablization transformation
            Y_n = sqrt(Y_n);
        end
        if  sum(Y_n)==0
            hyperparam_sigma_n = 1;
        else
            hyperparam_sigma_n = std(Y_n);
        end
        
        hyperparam_p_connected = output_warm(j).alpha;
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
        %output(j).w_estimate = alpha_tmp.*mu_tmp;
    end
    %toc
    %-------------------------------------------%
    time_record(2*t+1) = toc;
end

%% Waste code:
% num_this_batch = 2*max_x_group;
% index_this_batch = 2*max_x_group*(t-1)+ (1:2*max_x_group);
% start_this_batch = 2*max_x_group*(t-1);
% for i = 1:nquantile  % permuting the columns
%     perm_index(:,i) = x_group{i}(randperm(max_x_group));
% end
% locations_trials(index_this_batch,:)=[perm_index(:, 2*(1:num_sources)-1); perm_index(:, 2*(1:num_sources))];
% X_next = zeros(size(pi_dense_local,1), num_this_batch);
% X_next_all = zeros(size(pi_dense_all,1), num_this_batch);
% for l = 1:num_this_batch
%     locations_next = locations_trials(start_this_batch+l,:);
%     X_next(:,l) = min(.95,sum(pi_Z_local(:, locations_next) ,2));
%     X_next_all(:,l) = min(.95,sum(pi_Z_all(:,locations_next),2));
% end
% X_g(index_this_batch,:)  = X_next';
