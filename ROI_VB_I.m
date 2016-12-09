%% ROI detection: Stage I
%  Directly Stimulating cell somas

%% Find the neurons that locate in the region:
% 
% all_neuron_locations record locations of the neurons and their layer info
all_neuron_locations =zeros(0,4);
for i = 1:num_layers
    all_neuron_locations = [all_neuron_locations; [neuron_locations{i}(:,1:3) i*ones(size(neuron_locations{i},1),1) ]];
end

epsilon = 0.01;
x_low = postsyn_position(1) - region_width/2-epsilon;
x_upp = postsyn_position(1) + region_width/2;
y_low = postsyn_position(2) - region_height/2-epsilon;
y_upp = postsyn_position(2) + region_height/2;

% Identify neurons within this location:
neuron_in_region = zeros(size(all_neuron_locations,1),1); 
for i = 1:size(all_neuron_locations,1)
    if all_neuron_locations(i,1) > x_low & all_neuron_locations(i,1) < x_upp
        if all_neuron_locations(i,2) > y_low & all_neuron_locations(i,2) < y_upp
            neuron_in_region(i)=1;
        end
    end
end

Z = zeros(sum(neuron_in_region),3);
count = 1;
for i = 1:size(all_neuron_locations,1)
    if neuron_in_region(i) > 0
        Z(count,:) = [all_neuron_locations(i,1:2) postsyn_position(3)];
        count = count + 1;
    end 
end


%% Define the combination of trials 
% We want to stimulate neurons that are far apart from each other

% To do this, we split the cells by the quantiles of their x and y
% coordinates, then alternating between x and y.

% Find the quantiles:
nquantile = num_sources*2;

[~, ~, x_ranks] = unique(Z(:,1));
x_freq = x_ranks / size(Z,1);
[~, ~, y_ranks] = unique(Z(:,2));
y_freq = y_ranks / size(Z,1);

x_index = ceil(x_freq / (1/nquantile));
y_index = ceil(y_freq / (1/nquantile));

x_group = cell(nquantile,1);
y_group = cell(nquantile,1);
for i = 1:size(Z,1)
    x_group{x_index(i)} = [x_group{x_index(i)} i];
    y_group{y_index(i)} = [y_group{y_index(i)} i];
end
max_x_group = 0;
max_y_group = 0;
for i = 1:nquantile
    max_x_group = max(max_x_group, size(x_group{i},2));
    max_y_group = max(max_y_group, size(y_group{i},2));
end

joint_n = 2*(max_x_group + max_y_group); 
M = ceil(N/joint_n);

trial_locations_on_grid = zeros(N, num_sources);
count = 0;
for m = 1:M
    % Alternating between rows and columns!
    
    perm_index = ones(max_x_group,nquantile);
    for i = 1:nquantile  % permuting the columns
        n_this_group = size(x_group{i},2);
        perm_index(1:n_this_group,i) = x_group{i}(randperm(size(x_group{i},2)));
    end
    trials_m = [perm_index(:, 2*(1:num_sources)-1); perm_index(:, 2*(1:num_sources))];
    
    perm_index = ones(max_y_group,nquantile);
    for i = 1:nquantile  % permuting the columns
        n_this_group = size(y_group{i},2);
        perm_index(1:n_this_group,i) =  y_group{i}(randperm(size(y_group{i},2)));
    end
    trials_m = [trials_m; [perm_index(:, 2*(1:num_sources)-1); perm_index(:, 2*(1:num_sources))]];
    
    
    if N < (m*joint_n) 
        trial_locations_on_grid( ((m-1)*joint_n+1 ):N,:) = trials_m(1: (N-((m-1)*joint_n )),:);
    else
        trial_locations_on_grid( (m-1)*joint_n + (1:joint_n),:) = trials_m;
    end 
end

%% Calculate the light-induced probability 
% these parameters govern the time delay, as a function of the
% point-spread-function stimuli for a particular trial
% in seconds
pi_k = zeros(N,size(all_locations,1));
B = diag(all_locations(:,1:3)*inv(A)*all_locations(:,1:3)')*ones(1,num_sources);

for n = 1:N
   this_trial = trial_locations_on_grid(n,:);
   this_trial_locations = Z(this_trial,:);
   diff_mat=B+(ones(size(all_locations,1),1)*diag(this_trial_locations*inv(A)*this_trial_locations')')-2*all_locations(:,1:3)*inv(A)*this_trial_locations';

   pi_kr = exp(-0.5*(diff_mat));
   pi_k(n,:) = min(.95,sum(pi_kr,2));
end

pi_k_spike = pi_k;
pi_k_spike(pi_k_spike > .65) = 1; 

% firing delay means and variances
d_mean_nk = d_mean0 + (1.5 - pi_k)*d_mean_coef;
d_sigma_nk = d_sigma0 + (1 - pi_k)*d_sigma_coef;

% sample "ground truth" firing delay
D = exprnd(d_mean_nk/data_params.dt) + evoked_params.stim_start;

% sample "ground truth" stimulations
X = rand(N,K) < pi_k_spike; %.2 
X(D > 2000) = 0;

%% Generate a response, given D, Pi, X, w

Y = zeros(N,data_params.T);

for n = 1:N
    firing_neurons = X(n,:) & all_amplitudes' > 0;
    
    if any(all_amplitudes(firing_neurons) > 0)
        
        evoked_params.times = D(n,firing_neurons);
        evoked_params.a = all_amplitudes(firing_neurons);
        evoked_params.tau_r = all_tau_rise(firing_neurons)/data_params.dt;
        evoked_params.tau_f = all_tau_fall(firing_neurons)/data_params.dt;
    else
        evoked_params.times = [];
        evoked_params.a = [];
        evoked_params.tau_r = [];
        evoked_params.tau_f = [];
    end
    
    [Y(n,:), mpp_n] = gen_trace_noise(data_params,bg_params,evoked_params);
    if n == 1
        mpp = mpp_n;
    else
        mpp(n) = mpp_n;
    end
end







%% Calculate summary statistics 

related_mpp = mpp;

% only consider a small time window for events of interest
for i = 1:size(trial_locations_on_grid,1)
    if size(mpp(i).event_times,2) > 0
        indices = mpp(i).event_times>evoked_params.stim_start  & mpp(i).event_times< (400+evoked_params.stim_start);
        related_mpp(i).amplitudes = mpp(i).amplitudes(indices);
        related_mpp(i).event_times = mpp(i).event_times(indices);
    end 
end

covariates = zeros(size(trial_locations_on_grid,1), size(Z,1));
for i = 1:N
	covariates(i, trial_locations_on_grid(i,:)) = 1;    
end

% With intercept, it is rank deficient. 
% covariates_intercept = [ones(N,1) covariates];
% rank(covariates)
% size(covariates)
%% Dividing the events by their amplitudes
% Now divide the events by quantiles of the amplitudes 
% We use overlapped regions to avoid separation due to 
amplitude_threshold = quantile([related_mpp.amplitudes], (1/num_threshold)*[0:num_threshold]);
amp_related_count_trials = ones(size(trial_locations_on_grid,1),num_threshold-1);
for j = 1:(num_threshold-1)
    for i = 1:size(amp_related_count_trials,1)
        amp_related_count_trials(i,j) = sum(related_mpp(i).amplitudes>amplitude_threshold(j) & related_mpp(i).amplitudes<(amplitude_threshold(j+2)+0.01));
    end
end

%% Fit a rough model on data from Stage I



K_z = size(Z,1);

% Prior mean of mus
hyperparams_eta = zeros(K_z+1,1);
% Prior sd of mus
hyperparams_sigma_s = 1*ones(K_z+1,1);

pi_kr = exp(-0.5*squareform(pdist(Z(:,1:3),'mahalanobis',A)).^2);

pi_nk = zeros(N,K_z);
for n = 1:N
    pi_nk(n,:) = min(1,sum(pi_kr(:,trial_locations_on_grid(n,:)),2)');
end

output= struct([]);
for j = 1:size(amp_related_count_trials,2)
    
    Y_n = amp_related_count_trials(:,j);
    hyperparams_sigma_n = std(Y_n);
    
    hyperparam_p_connected = .1*ones(K_z+1,1);
    
    options= struct();
    options.verbose= true;
    options.center= 0;
    
    X_temp = pi_nk;
    X_temp = [ones(size(X_temp,1),1) X_temp];
    [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs_general(X_temp, Y_n, ...
        hyperparams_sigma_n, hyperparams_sigma_s, hyperparam_p_connected,...
        hyperparams_eta,options);
    
    output(j).alpha = alpha_tmp;
    output(j).mu = mu_tmp;
    output(j).s_sq = s_sq_tmp;
    
    output(j).w_estimate = alpha_tmp.*mu_tmp;
end

%------------------------End of first stage-------------------------------------%

%
