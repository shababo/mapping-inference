%% Option II:
% Fit a "full" model for evaluation
%%
% initialize alpha, mu, and s_sq to the estimates from the reduced model
fprintf('Evaluate Batch %d of %d\n', t, N);
%tic
% Designing trials
    end_this_batch = num_trials_first + (t-1)*num_trials_batch;

amplitude_threshold=  quantile([mpp_new.amplitudes], (1/num_threshold)*[0:num_threshold]);
% Count the events in each amplitude bins:
related_mpp_n=struct();
for l = 1:end_this_batch
    if size(mpp_new(l).event_times,2) > 0
        indices = mpp_new(l).event_times>evoked_params.stim_start ...
            & mpp_new(l).event_times< (400+evoked_params.stim_start);
        related_mpp_n(l).amplitudes = mpp_new(l).amplitudes(indices);
        related_mpp_n(l).event_times = mpp_new(l).event_times(indices);
    else
        related_mpp_n(l).amplitudes = [];
        related_mpp_n(l).event_times = [];
    end
    for j = 1:num_threshold
        Y_g(l,j) = sum(related_mpp_n(l).amplitudes>amplitude_threshold(j) & related_mpp_n(l).amplitudes<(amplitude_threshold(j+1)+0.01));
    end
end
%toc
%---------------------------------------------------------%
% Part III:
% Fit the VB model to update the parameters:
% Initialization:
output_warm = output_refit;
X_temp = X_g(1:end_this_batch,:);
X_temp = [ones(size(X_temp,1),1) X_temp];
for j = 1:size(Y_g,2)
    fprintf('%d',j);
    Y_n = Y_g(1:end_this_batch,j);
    % Variance stablization transformation
   if sqrt_transform 
    Y_n = sqrt(Y_n);
   end
   
    if  sum(Y_n)==0
        hyperparam_sigma_n = 1;
    else
        hyperparam_sigma_n = std(Y_n);
    end
    
	
    hyperparam_p_connected = .1*ones(K_z+1,1);
    hyperparam_eta =  zeros(K_z+1,1);
    hyperparam_sigma_s =  ones(K_z+1,1);
    
    options = struct();
    options.alpha =  output_warm(j).alpha;
    options.mu= output_warm(j).mu;
    options.verbose= false;
    options.center= 0;
    [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs_general(X_temp, Y_n, ...
        hyperparam_sigma_n, hyperparam_sigma_s, hyperparam_p_connected, ...
        hyperparam_eta, options);
    output_refit(j).alpha = alpha_tmp;
    output_refit(j).mu = mu_tmp;
    output_refit(j).s_sq = s_sq_tmp;
    %output(j).w_estimate = alpha_tmp.*mu_tmp;
end
%toc
%-------------------------------------------%

