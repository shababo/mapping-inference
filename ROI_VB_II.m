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
sigma_s = 1;
sigma_n=1;
num_samples=50; % Number of samples to estimate the expected entropies 
num_trials_batch=20;
%% Design parameters:

ND=100; % Number of batches

num_simultaneous_inputs = 1; % Number of locations to stimulate in each trial

% We consider new locations on a grid (rather than the cell somas) in order
% to gain extra information (e.g., to be able to distinguish two
% neighbouring cells)

num_dense=50; % Grid density 

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
    size(selected_cells{j})
    alphas = zeros(size(selected_cells{j},1),1); %ones(params.K, 1) * alpha_0;
    mu = zeros(size(selected_cells{j},1), 1);
    s_sq = zeros(size(selected_cells{j},1),1);
    options= struct();
    options.verbose= false;
    options.center= 0;
    
    for sample = 1:n_varbvs_samples
        % Prepare the data matrix:
        X_temp = pi_nk_j>rand(params.N,size(pi_nk_j,2));
        [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs_general(X_temp, Y_n, hyperparam_sigma_n, params.sigma_s(1), hyperparam_p_connected(1), params.eta,options);
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
output_warm= output_reduced;

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
tic
locations_record = [];    
for t = 1:ND
    
    fprintf('trial %d of %d\n', t, ND);
    %--------------------------------------------------%
    % Part I: evaluate the entropies
    
    delta_H = zeros(num_dense^2, size(amp_related_count_trials,2));
    for j = 1:size(amp_related_count_trials,2)
        % compute current objective function value
        H_current = per_neuron_joint_entropy(output_greedy(j).alpha, output_greedy(j).mu, output_greedy(j).s_sq);
        % for each neuron, calculate the expected change in the objective function
        H_expected = approximate_expected_joint_entropy_single_neuron(output_greedy(j).alpha, output_greedy(j).mu, ...
            output_greedy(j).s_sq, hyperparam_sigma_n, num_samples);
        delta_H_cell = H_current-H_expected;  
        % Map the change of entropies to the dense grid 
        % ** Switch to a search algorithm. 
        pi_kr_j = pi_kr_grid(output_reduced(j).cell_list,output_reduced(j).locations);
        delta_H(output_reduced(j).locations,j) = pi_kr_j' * delta_H_cell;  
    end
    
    % Display the heatmap:
    heat_all = sum(delta_H,2);
    % reformat into matrix
    heat_mat = zeros(num_dense,num_dense);
    
    for i = 1:num_dense
        for l = 1:num_dense
            heat_mat(i,l)= heat_all((i-1)*num_dense + l,:);
        end
    end
    [cent, varargout]=FastPeakFind(heat_mat,0);
   
%     HeatMap(heat_mat)
%       xlim([0,50]);
%     ylim([0,50]);
% 
% 
%     HeatMap(varargout)
%   xlim([0,50]);
%     ylim([0,50]);

    peaks = zeros(size(cent,1)/2,2);
    for i = 1:(size(cent,1)/2)
        peaks(i,:) = cent([1+ 2*(i-1) 2*(i-1)+2]);
    end
    
    idx= (peaks(:,2)-1)*num_dense + peaks(:,1);
    %---------------------------------------------%
    
    % Part II: choose the locations and generate new data 
    %[~, ix] = sort(sum(delta_H,2), 1, 'descend');
    
    % Find multiple treaments: 
    X_next=[];
    X_next_all=[];
    for l = 1:num_trials_batch
        locations_next = randsample(size(idx,1), num_simultaneous_inputs);
        locations_record = [locations_record idx(locations_next)];
        X_next = [X_next, sum(pi_kr_grid(:, idx(locations_next)) ,2)];
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
    parfor j = 1:size(Y_g,2)
        n_selected_cells=size(output_greedy(j).cell_list,2);
        Y_n = Y_g(:,j);
        hyperparam_sigma_n = sqrt(length(params.t))*params.sigma_n/abs(alpha_sum);
        hyperparam_p_connected = .1*ones(n_selected_cells,1);
        
        alphas = zeros(n_selected_cells,1); %ones(params.K, 1) * alpha_0;
        mu = zeros(n_selected_cells, 1);
        s_sq = zeros(n_selected_cells,1);
        options = struct();
        options.alpha =  output_warm(j).alpha;
        options.mu= output_warm(j).mu;
        options.verbose= false;
       options.center= 0;
      
       % Drop the sampling step, run regression on the probability instead
       
       n_varbvs_samples = 1;
        %for sample = 1:n_varbvs_samples
            %X_temp = X_g(:,output_greedy(j).cell_list)>rand(params.N+t*num_trials_batch, n_selected_cells);
            X_temp = X_g(:,output_greedy(j).cell_list);
            [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs_general(X_temp, Y_n, hyperparam_sigma_n, params.sigma_s(1), hyperparam_p_connected(1), params.eta, options);
            alphas = alphas+alpha_tmp/n_varbvs_samples;
            mu = mu+mu_tmp/n_varbvs_samples;
            s_sq = s_sq+s_sq_tmp/n_varbvs_samples;
        %end
        
        output_greedy(j).alpha = alphas;
        output_greedy(j).mu = mu;
        output_greedy(j).s_sq = s_sq;
        output_greedy(j).w_estimate = alphas.*mu;
        
    end
    output_warm = output_greedy;
    %-------------------------------------------%
    
    
    % calculate metrics
    %SP_greedy(t) = sum(gamma ~= (alpha_greedy >= .5));
    %log_loss_greedy(t) = log_loss(gamma, alpha_greedy);
    %rmse_greedy(t) = rms_error(alpha_greedy, mu_greedy, w);
    
    % save
    %         if mod(t, 100) == 0
    %     if t == ceil(T/2)
    %         disp('saving temporary output to disk');
    %         save(sprintf(outfile_path, t));
    %     end
    
end
toc


%% Fit a semi-full model to summarize the result
output_warm = output;
output_post= struct([]);
tic
for j = 1:size(amp_related_count_trials,2)
    
     Y_n = Y_g(:,j);
       
    hyperparam_sigma_n = sqrt(length(params.t))*params.sigma_n/abs(alpha_sum);
    
    hyperparam_p_connected = .1*ones(params.K,1);
    
    alphas = zeros(params.K,1); %ones(params.K, 1) * alpha_0;
    mu = zeros(params.K, 1);
    s_sq = zeros(params.K,1);
    n_varbvs_samples = 5;
      options = struct();
        options.alpha =  output_warm(j).alpha;
        options.mu= output_warm(j).mu;
        options.verbose= false;
       options.center= 0;
      
    for sample = 1:n_varbvs_samples
        X_temp = X_g>rand(size(X_g,1), size(X_g,2));
        [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs_general(X_temp, Y_n, hyperparam_sigma_n, params.sigma_s(1), hyperparam_p_connected(1), params.eta,options);
        alphas = alphas+alpha_tmp/n_varbvs_samples;
        mu = mu+mu_tmp/n_varbvs_samples;
        s_sq = s_sq+s_sq_tmp/n_varbvs_samples;
    end
    
    output_post(j).alpha = alphas;
    output_post(j).mu = mu;
    output_post(j).s_sq = s_sq;
    output_post(j).w_estimate = alphas.*mu;
end
toc

%% Visualization:

figure(80)

for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*25);
    set(temp,'MarkerFaceColor','k');
   alpha(temp,0.3);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

%selected_pixels = zeros(num_dense,num_dense,size(amp_related_count_trials,2));
%cent = cell(size(amp_related_count_trials,2),1);

    xlim([20,460]);
    ylim([-900,-400]);
%
%      potential_neuron_grid = scatter(Z(:,1),...
%      -Z(:,2),20,colormap(2,:),...
%     'filled','d');
%     set(potential_neuron_grid,'MarkerFaceColor','k');
%     alpha(potential_neuron_grid,0.2);

for j = 1:size(amp_related_count_trials,2)
    coef = output(j).alpha;
    coef_thres = 0.05;% quantile(coef,0.98);
    potential_neuron_grid = scatter(Z(coef>coef_thres,1), -Z(coef>coef_thres,2), amplitude_threshold(j+1)*25, 'filled','o');
    set(potential_neuron_grid,'MarkerFaceColor','r');
    alpha(potential_neuron_grid,0.3);
    hold on

end

hold off
view(2)

%%

figure(100)

for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*25);
    set(temp,'MarkerFaceColor','k');
   alpha(temp,0.3);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

%selected_pixels = zeros(num_dense,num_dense,size(amp_related_count_trials,2));
%cent = cell(size(amp_related_count_trials,2),1);

    xlim([20,460]);
    ylim([-900,-400]);


for j = 1:size(amp_related_count_trials,2)
    coef = output_post(j).alpha;
    coef_thres = 0.05;
    potential_neuron_grid = scatter(Z(coef>coef_thres,1), ...
        -Z(coef>coef_thres,2), amplitude_threshold(j+1)*25,'filled','o');
    set(potential_neuron_grid,'MarkerFaceColor','r');
    alpha(potential_neuron_grid,0.4);
    hold on

end

    stimuli_grid = scatter(Z_dense(locations_record,1), ...
        -Z_dense(locations_record,2),10);
    set(stimuli_grid,'MarkerFaceColor','g');
    alpha(stimuli_grid,1);


hold off
view(2)
