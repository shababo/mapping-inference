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
num_samples=50; % Number of samples to estimate the expected entropies
num_trials_batch=20;
%% Design parameters:

ND=100; % Number of batches

num_simultaneous_inputs = 4; % Number of locations to stimulate in each trial
num_peaks = 30;
% We consider new locations on a grid (rather than the cell somas) in order
% to gain extra information (e.g., to be able to distinguish two
% neighbouring cells)

num_dense=80; % Grid density

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

%% Option I: evaluate the entropy on a reduced model
% VB inference in each batch is very fast
% Its validity relies on sufficient data from the first stage 


% Given the VB estimates, we can reduce the number of cells to consider for
% each amplitude bin by removing those with probability of firing being
% in the lower 80% (purely to reduce the dimension)

% The reduced models have far fewer parameters to fit.

output_reduced= struct([]);
selected_cells = cell(size(amp_related_count_trials,2),1);
idx_temp = 1:params.K;
for j = 1:size(amp_related_count_trials,2)
    
    Y_n = amp_related_count_trials(:,j);
    params.sigma_n = std(Y_n);
    % Select the cells whose probability of excitation is larger than 0.05
    selected_cells{j} = idx_temp(output(j).alpha > 0.1);
    pi_nk_j =pi_nk(:,selected_cells{j});
    %size(selected_cells{j})
    hyper_sigma_s = 1*ones( size(selected_cells{j},2)+1,1);
    hyperparam_p_connected = .1*ones( size(selected_cells{j},2)+1,1);
    
    alphas = zeros(size(selected_cells{j},2)+1,1); %ones(params.K, 1) * alpha_0;
    mu = zeros(size(selected_cells{j},2)+1, 1);
    s_sq = zeros(size(selected_cells{j},2)+1,1);
    options= struct();
    options.verbose= false;
    options.center= 0;
    
    X_temp = pi_nk_j;
    X_temp = [ones(size(X_temp,1),1) X_temp];
    [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs_general(X_temp, Y_n, params.sigma_n, hyper_sigma_s, hyperparam_p_connected, params.eta,options);
    output_reduced(j).alpha = alpha_tmp;
    output_reduced(j).mu = mu_tmp;
    output_reduced(j).s_sq = s_sq_tmp;
    
    output_reduced(j).w_estimate = alpha_tmp.*mu_tmp;
    output_reduced(j).cell_list = selected_cells{j};
end
output_warm= output_reduced;

%% Calculate and save the innner products of weight matrix w
% The responsive cells might not cover the full space
% Hence, we can further save some computation by considering only the
% locations that are relevant (with total firing probability > 0.1)
cell_list = unique([output_reduced.cell_list]);
weights_mat_full = pi_kr_grid(cell_list, :);

% Calculate the inner product of the induced probability
inner_products = weights_mat_full'*weights_mat_full;

self_products = diag(inner_products)*ones(1,size(inner_products,1));
inner_normalized_products = inner_products./self_products;
%
% a= weights_mat_full(:,1);
% b= weights_mat_full(:,2);
%
% k= a'*b/b'*b
% a'*a /b'*b + k^2  - 2k a'*b/b'*b
% (a- k*b)'*(a- k*b)
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
tic
locations_record = [];
time_record = zeros(ND,1);
for t = 1:ND
    fprintf('trial %d of %d\n', t, ND);
    %--------------------------------------------------%
    % Part I: evaluate the entropies
     delta_H = zeros(params.K,1);
    for j = 1:size(amp_related_count_trials,2)
        % compute current objective function value
        H_current = per_neuron_joint_entropy(output_greedy(j).alpha(2:end), output_greedy(j).mu(2:end), output_greedy(j).s_sq(2:end));
        % for each neuron, calculate the expected change in the objective function
        hyperparam_sigma_n = std(amp_related_count_trials(:,j));
        
        H_expected = approximate_expected_joint_entropy_single_neuron(output_greedy(j).alpha(2:end), output_greedy(j).mu(2:end), ...
            output_greedy(j).s_sq(2:end), hyperparam_sigma_n, num_samples);
        delta_H(output_reduced(j).cell_list) = delta_H(output_reduced(j).cell_list) + H_current-H_expected;
    end
    entropy_locations = weights_mat_full'*delta_H(cell_list);
    
    idx = [];
    for i = 1:num_peaks
        [m_id temp_idx] = max(entropy_locations);
        idx = [idx temp_idx];
        entropy_locations = entropy_locations - inner_normalized_products(:,temp_idx)*m_id;
    end
    %---------------------------------------------%
    
    % Part II: choose the locations and generate new data
    %[~, ix] = sort(sum(delta_H,2), 1, 'descend');
    
    % Find multiple treaments:
    X_next=[];
    X_next_all=[];
    for l = 1:num_trials_batch
        locations_next = randsample(size(idx,2), num_simultaneous_inputs);
        % Add some
        locations_record = [locations_record idx(locations_next)];
        X_next = [X_next, min(.95,sum(pi_kr_grid(:, idx(locations_next)) ,2))];
        X_next_all = [X_next_all, min(.95,sum(pi_kr_all(:,idx(locations_next)),2))];
    end
    
    % update X matrix
    X_g = [X_g; X_next'];
    
    % Draw a new sample using the same mechanism in
    % gendata_fullmodel_multicell.m
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
        
        %Y_next = zeros(1,data_params.T);
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
        
        [Y_new( (t-1)*num_trials_batch+l,:), mpp_n] = gen_trace_noise(data_params,bg_params,evoked_params);
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
            amp_related_count_trials_new((t-1)*num_trials_batch+l,j) = sum(related_mpp_n(l).amplitudes>amplitude_threshold(j) & related_mpp_n(l).amplitudes<(amplitude_threshold(j+2)+0.01));
        end
        
    end
    
    Y_g = [Y_g; amp_related_count_trials_new( (t-1)*num_trials_batch+(1:num_trials_batch),:) ];
    
    %---------------------------------------------------------%
    % Part III:
    % Fit the VB model to update the parameters:
    % Note: if there is no events at all, the algorithm will report a
    % warning
    output_warm = output_greedy;
    fprintf('Sum of sq %d \n', sum(output_warm(1).s_sq));
    
    for j = 1:size(Y_g,2)
       % fprintf('%d',j);
        n_selected_cells=size(output_greedy(j).cell_list,2);
        Y_n = Y_g((t-1)*num_trials_batch+(1:num_trials_batch),j);
        if  sum(Y_n)==0 
        hyperparam_sigma_n = 1;
        else
        hyperparam_sigma_n = std(Y_n);
            
        end
        
        hyperparam_p_connected = [output_warm(j).alpha];
        
        % Prior distribution for mu:
        hyperparam_eta =  [output_warm(j).mu];
        hyperparam_sigma_s = sqrt(output_warm(j).s_sq);
        
        % fprintf('Sum of mu %d \n', sum(output_warm(j).mu));
        alphas = zeros(n_selected_cells+1,1); %ones(params.K, 1) * alpha_0;
        mu = zeros(n_selected_cells+1, 1);
        s_sq = zeros(n_selected_cells+1,1);
        options = struct();
        options.alpha =  output_warm(j).alpha;
        options.mu= output_warm(j).mu;
        options.verbose= false;
        options.center= 0;
        
        
        % Drop the sampling step, run regression on the probability instead
        %for sample = 1:n_varbvs_samples
        %X_temp = X_g(:,output_greedy(j).cell_list)>rand(params.N+t*num_trials_batch, n_selected_cells);
        X_temp = X_g((t-1)*num_trials_batch+(1:num_trials_batch),output_greedy(j).cell_list);
        X_temp = [ones(size(X_temp,1),1) X_temp];
        
        [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs_general(X_temp, Y_n, hyperparam_sigma_n, hyperparam_sigma_s, hyperparam_p_connected, hyperparam_eta, options);
        alphas = alphas+alpha_tmp/n_varbvs_samples;
        mu = mu+mu_tmp/n_varbvs_samples;
        s_sq = s_sq+s_sq_tmp/n_varbvs_samples;
        %end
        
        output_greedy(j).alpha = alphas;
        output_greedy(j).mu = mu;
        output_greedy(j).s_sq = s_sq;
        output_greedy(j).w_estimate = alphas.*mu;
        
    end
    
    %-------------------------------------------%
    time_record(t) = toc;
end
