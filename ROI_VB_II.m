%% Loading functions and Data generation
clear;
addpath(genpath('../psc-detection'),genpath('../mapping-inference'),genpath('../mapping-core'));
%addpath(genpath('../Data'));
%%
% Generate data 
run('gendata_fullmodel_multicells.m')
% Run Stage I 
run('ROI_VB_I.m')


%% Hyperparameters:
% Might b
alpha = 0.1;
sigma_s = 1;
sigma_n=1;
num_samples=5; % Number of samples to estimate the expected entropies 
%% Design parameters:

ND=1; % Number of trials to design

num_simultaneous_inputs = 4; % Number of locations to stimulate in each trial

% We consider new locations on a grid (rather than the cell somas) in order
% to gain extra information (e.g., to be able to distinguish two
% neighbouring cells)

num_dense=20; % Grid density 

% The new stimuli are chosen to stay away from the true neurons 
% We stimulate a dense grid instead
x_dense = zeros(num_dense,1);
y_dense = zeros(num_dense,1);
x_dense = (0:(num_dense-1))*(max(Z(:,1))-min(Z(:,1)))/(num_dense-1) + min(Z(:,1));
y_dense = (0:(num_dense-1))*(max(Z(:,2))-min(Z(:,2)))/(num_dense-1) + min(Z(:,2));
Z_dense = zeros(num_dense^2,3);

 for i = 1:num_dense
        for l = 1:num_dense
            Z_dense((i-1)*num_dense + l,:) = [x_dense(i) y_dense(l) postsyn_position(3)];
        end
 end

% Calculate the probability of firing for ALL neurons 
% We will use it in simulating the real spikes 
B = diag(all_locations*inv(A)*all_locations');
diff_mat=B*ones(size(Z_dense,1),1)'+(ones(size(all_locations,1),1)*diag(Z_dense*inv(A)*Z_dense')')-2*all_locations*inv(A)*Z_dense';
pi_kr_all = exp(-0.5*(diff_mat));


% Calculate the probability of firing for the LOCAL neurons (i.e., those
% that are within the 2-D plane we consider)
% We will use this in estimating the expectation, and in fitting the model 
B = diag(params.coords*inv(A)*params.coords');
diff_mat=B*ones(size(Z_dense,1),1)'+(ones(size(params.coords,1),1)*diag(Z_dense*inv(A)*Z_dense')')-2*params.coords*inv(A)*Z_dense';
pi_kr_grid = exp(-0.5*(diff_mat));

%% Fit a reduced model 
% Given the VB estimates, we can reduce the number of cells to consider for
% each amplitude bin by removing those with probability of firing being
% less than 0.05 (or any threshold)

% The reduced models have far fewer parameters to fit. 

output_reduced= struct([]);
selected_cells = cell(size(amp_related_count_trials,2),1);
idx_temp = 1:params.K;
for j = 1:size(amp_related_count_trials,2)
    
    Y_n = amp_related_count_trials(:,j);
    hyperparam_sigma_n = sqrt(length(params.t))*params.sigma_n/abs(alpha_sum);
    hyperparam_p_connected = .1*ones(params.K,1);
    
    n_varbvs_samples = 5;

    % Select the cells whose probability of excitation is larger than 0.05
    selected_cells{j} = idx_temp(output(j).alpha >0.05);
    pi_nk_j =pi_nk(:,selected_cells{j});
    
    alphas = zeros(size(selected_cells{j},1),1); %ones(params.K, 1) * alpha_0;
    mu = zeros(size(selected_cells{j},1), 1);
    s_sq = zeros(size(selected_cells{j},1),1);
    
    for sample = 1:n_varbvs_samples
        % Prepare the data matrix: 
        X_temp = pi_nk_j>rand(params.N,size(pi_nk_j,2));
   
        [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs_general(X_temp, Y_n, hyperparam_sigma_n, params.sigma_s(1), hyperparam_p_connected(1), params.eta);
        alphas = alphas+alpha_tmp/n_varbvs_samples;
        mu = mu+mu_tmp/n_varbvs_samples;
        s_sq = s_sq+s_sq_tmp/n_varbvs_samples;
    end
    
    output_reduced(j).alpha = alphas;
    output_reduced(j).mu = mu;
    output_reduced(j).s_sq = s_sq;
    
    output_reduced(j).w_estimate = alphas.*mu;
    output_reduced(j).cell_list = selected_cells{j};
end

%% Selected the relevant locations for each amplitude bins
% The responsive cells might not cover the full space
% Hence, we can further save some computation by considering only the
% locations that are relevant (with total firing probability > 0.1)

idx_temp = 1:size(Z_dense,1);
for j = 1:size(amp_related_count_trials,2)
    pi_kr_grid_j = pi_kr_grid(output_reduced(j).cell_list,:); 
    output_reduced(j).locations = idx_temp( sum(pi_kr_grid_j,1)>0.1);
end
%% Select the optimal location(s)

% It takes ~70 seconds for each t..

% Storage:
H_actual_history = NaN(ND,num_threshold);
H_expected_next_history = NaN(num_threshold,num_dense^2);


% Initialize the response and covariate using the batched data 
X_g = pi_nk;
Y_g = amp_related_count_trials;

% Storage for new data (might be unnecessary)
trace_new = zeros(ND, data_params.T);
amp_related_count_trials_new = ones(ND,num_threshold-1);

% initialize alpha, mu, and s_sq to the estimates from the reduced model 
output_greedy = output_reduced;
obj_function = @joint_sparsity_weight_entropy; % 

for t = 1:ND
    
    fprintf('trial %d of %d\n', t, ND);
    tic
    %--------------------------------------------------%
    % Part I: evaluate the entropies
    
    H_expected = zeros(num_dense^2,1);
    delta_H = zeros(num_dense^2, size(amp_related_count_trials,2));
    
    for j = 1:size(amp_related_count_trials,2)
        
        % compute current objective function value
        H_current = obj_function(output_greedy(j).alpha, output_greedy(j).mu, output_greedy(j).s_sq);
        
        % for each neuron, calculate the expected change in the objective function
        X_g_j = X_g(:, output_greedy(j).cell_list);
        for k = output_reduced(j).locations
            X_t_plus_1 = [X_g_j; zeros(1,size(output_greedy(j).cell_list,2))];
            X_t_plus_1(end, :) = pi_kr_grid(output_reduced(j).cell_list,k);
            H_expected(k) = compute_expectation_fast2(num_samples, obj_function, output_greedy(j).alpha, ...
                output_greedy(j).mu, output_greedy(j).s_sq, sigma_n, sigma_s, alpha, X_t_plus_1, Y_g(:,j),params.eta);
            X_t_plus_1(end, k) = 0;
            disp([num2str(k) ' evaluated.']);
        end
        delta_H(:,j) = H_current - H_expected;
    end
    %---------------------------------------------%
    
    % Part II: choose the locations and generate new data 
    [~, ix] = sort(sum(delta_H,2), 1, 'descend');
    X_next = sum(pi_kr_grid(:,ix(1:num_simultaneous_inputs)),2);
    
    % update X matrix
    X_g = [X_g; X_next'];
    
    % Draw a new sample using the same mechanism in
    % gendata_fullmodel_multicell.m 
    X_next_all = sum(pi_kr_all(:,ix(1:num_simultaneous_inputs)),2);
    pi_k_spike = X_next_all;
    pi_k_spike(pi_k_spike > .65) = 1;
    
    % firing delay means and variances
    d_mean_nk = d_mean0 + (1.5 - X_next_all)*d_mean_coef;
    d_sigma_nk = d_sigma0 + (1 - X_next_all)*d_sigma_coef;
    
    % sample "ground truth" firing delay
    D = exprnd(d_mean_nk/data_params.dt) + evoked_params.stim_start;
    
    % sample "ground truth" stimulations
    X = rand(size(D,1),1) < pi_k_spike; %.2
    X(D > 2000) = 0;
    
    Y_next = zeros(1,data_params.T);
    firing_neurons = X & all_amplitudes > 0;
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
    
    [Y_new(t,:), mpp_n] = gen_trace_noise(data_params,bg_params,evoked_params);
    if t == 1
        mpp_new = mpp_n;
    else
        mpp_new(t) = mpp_n;
    end
    
    
    % Count the events in each amplitude bins:
    related_mpp_n=mpp_n;
    if size(mpp_n.event_times,2) > 0
        indices = mpp_n.event_times>evoked_params.stim_start  & mpp_n.event_times< (400+evoked_params.stim_start);
        related_mpp_n.amplitudes = mpp_n.amplitudes(indices);
        related_mpp_n.event_times = mpp_n.event_times(indices);
    end
    
    for j = 1:(num_threshold-1)
        amp_related_count_trials_new(t,j) = sum(related_mpp_n.amplitudes>amplitude_threshold(j) & related_mpp_n.amplitudes<(amplitude_threshold(j+2)+0.01));
    end
    
    Y_g = [Y_g; amp_related_count_trials_new(t,:)];
    
    %---------------------------------------------------------%
    
    % Part III:
    % Fit the VB model to update the parameters:
    for j = 1:size(Y_g,2)
        n_selected_cells=size(output_greedy(j).cell_list,1);
        Y_n = Y_g(:,j);
        hyperparam_sigma_n = sqrt(length(params.t))*params.sigma_n/abs(alpha_sum);
        hyperparam_p_connected = .1*ones(n_selected_cells,1);
        
        n_varbvs_samples = 5;
        alphas = zeros(n_selected_cells,1); %ones(params.K, 1) * alpha_0;
        mu = zeros(n_selected_cells, 1);
        s_sq = zeros(n_selected_cells,1);
        for sample = 1:n_varbvs_samples
            X_temp = X_g(:,output_greedy(j).cell_list)>rand(params.N+t,n_selected_cells);
            [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs_general(X_temp, Y_n, hyperparam_sigma_n, params.sigma_s(1), hyperparam_p_connected(1), params.eta);
            alphas = alphas+alpha_tmp/n_varbvs_samples;
            mu = mu+mu_tmp/n_varbvs_samples;
            s_sq = s_sq+s_sq_tmp/n_varbvs_samples;
        end
        
        output_reduced(j).alpha = alphas;
        output_reduced(j).mu = mu;
        output_reduced(j).s_sq = s_sq;
        
        output_reduced(j).w_estimate = alphas.*mu;
    end
    %-------------------------------------------%
    
    
    % calculate metrics
    %SP_greedy(t) = sum(gamma ~= (alpha_greedy >= .5));
    %log_loss_greedy(t) = log_loss(gamma, alpha_greedy);
    %rmse_greedy(t) = rms_error(alpha_greedy, mu_greedy, w);
    
    toc
    % save
    %         if mod(t, 100) == 0
    %     if t == ceil(T/2)
    %         disp('saving temporary output to disk');
    %         save(sprintf(outfile_path, t));
    %     end
    
end


