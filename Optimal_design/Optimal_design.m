%% Loading functions and Data generation
clear;
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Parameters for the data generating mechanism
rng(12242,'twister');
%rng(12212,'twister');

% load parameters for the model
run('./Parameters_setup_ground_truth.m')

%% 
% Pre-processing for  experiment
num_samples=50; % Number of samples to estimate the expected entropies
num_trials_batch=20;

num_peaks = 20;
% We consider new locations on a grid (rather than the cell somas) in order
% to gain extra information (e.g., to be able to distinguish two
% neighbouring cells)

num_dense=40; % Grid density

N=200; % Number of batches
%N_initial = 20;
num_trials_first = 200;

%forget_scale = 0.5;
sqrt_transform = false;
%design = 'random'; 
num_threshold=10;
%%
% covariance of point spread function
A = diag([200, 200, 750]);
for randomseed = 1:10    
    flnm=strcat('../Data/sim-results/A200S', num2str(randomseed),'.mat'); 
    num_sources = 4;  % Number of locations to stimulate in each trial
    run('Parameters_setup_experiment.m')
    % Run analysis and design
    
    design = 0;
    output= struct([]);
    for j = 1:num_threshold
        output(j).alpha = .1*ones(K_z+1,1);
        output(j).mu = zeros(K_z+1,1);
        output(j).s_sq = ones(K_z+1,1);
        %output(j).w_estimate = zeros(K_z+1,1);
    end
    output_refit = output;
    obj_function = @joint_sparsity_weight_entropy; %
    X_g = zeros((N-1)*num_trials_batch + num_trials_first,K_z);
    locations_trials = zeros( (N-1)*num_trials_batch + num_trials_first,num_sources);
    Y_g = zeros((N-1)*num_trials_batch + num_trials_first,num_threshold);
    grid_index = 1:size(pi_dense_all,2);
    
    tic
    locations_record = [];
    time_record = zeros(2*N+1,1);
    time_record(1) = toc;
    
    output_random =cell(N/10,1);
    output_random_refit =cell(N/10,1);
    for outer_i = 1:(N/10)
        t_start = (outer_i-1)*10+1;
        t_end = outer_i*10;
        run('Experiment_full.m')
        output_random{outer_i} = output;
        run('Experiment_refit.m')
        output_random_refit{outer_i} = output_refit;
    end
    X_random = X_g;
    mpp_random = mpp_new;
    time_record_random = time_record;
    location_random = locations_trials;
    
    % Optimal design
    
    design = 1;
    output= struct([]);
    for j = 1:num_threshold
        output(j).alpha = .1*ones(K_z+1,1);
        output(j).mu = zeros(K_z+1,1);
        output(j).s_sq = ones(K_z+1,1);
        %output(j).w_estimate = zeros(K_z+1,1);
    end
    output_refit = output;
    obj_function = @joint_sparsity_weight_entropy; %
    X_g = zeros((N-1)*num_trials_batch + num_trials_first,K_z);
    locations_trials = zeros((N-1)*num_trials_batch + num_trials_first,num_sources);
    Y_g = zeros((N-1)*num_trials_batch + num_trials_first,num_threshold);
    grid_index = 1:size(pi_dense_all,2);
    
    tic
    locations_record = [];
    time_record = zeros(2*N+1,1);
    time_record(1) = toc;
    
    output_optimal =cell(N/10,1);
    output_optimal_refit =cell(N/10,1);
    for outer_i = 1:(N/10)
        t_start = (outer_i-1)*10+1;
        t_end = outer_i*10;
        run('Experiment_full.m')
        output_optimal{outer_i} = output;
        run('Experiment_refit.m')
        output_optimal_refit{outer_i} = output_refit;
    end
    X_optimal = X_g;
    mpp_optimal = mpp_new;
    time_record_optimal = time_record;
    location_optimal = locations_trials;
    
    % near othogonal design
    sigmoid_a=1;
    design = 3;
    output= struct([]);
    for j = 1:num_threshold
        output(j).alpha = .1*ones(K_z+1,1);
        output(j).mu = zeros(K_z+1,1);
        output(j).s_sq = ones(K_z+1,1);
        %output(j).w_estimate = zeros(K_z+1,1);
    end
    output_refit = output;
    obj_function = @joint_sparsity_weight_entropy; %
    X_g = zeros((N-1)*num_trials_batch + num_trials_first,K_z);
    locations_trials = zeros((N-1)*num_trials_batch + num_trials_first,num_sources);
    Y_g = zeros((N-1)*num_trials_batch + num_trials_first,num_threshold);
    grid_index = 1:size(pi_dense_all,2);
    
    tic
    locations_record = [];
    time_record = zeros(2*N+1,1);
    time_record(1) = toc;
    
    output_ortho =cell(N/10,1);
    output_ortho_refit =cell(N/10,1);
    for outer_i = 1:(N/10)
        t_start = (outer_i-1)*10+1;
        t_end = outer_i*10;
        run('Experiment_full.m')
        output_ortho{outer_i} = output;
        run('Experiment_refit.m')
        output_optimal_refit{outer_i} = output_refit;
    end
    X_ortho = X_g;
    mpp_ortho = mpp_new;
    time_record_ortho = time_record;
    location_ortho = locations_trials;
    
    
    % Evaluate the performance: 
% 1, normalized reconstruction error of connectivity
% 2, AUC for connectivity
% 3, normalized amplitudes reconstruction
bandwidth = [200 200];
output_eva = output_random; 
run('Experiment_evaluate.m')

prob_inclusion_random = overall_connectivity;
weighted_inclusion_random = overall_connectivity_w;

NRE_conn_random = NRE_conn;
NRE_amp_random = NRE_amp;
AUC_conn_random = AUC_conn;
Spatial_amp_random = Spatial_amp;
Spatial_conn_random = Spatial_conn;
    

output_eva = output_optimal; 
run('Experiment_evaluate.m')

prob_inclusion_optimal = overall_connectivity;
weighted_inclusion_optimal= overall_connectivity_w;
NRE_conn_optimal = NRE_conn;
NRE_amp_optimal = NRE_amp;
AUC_conn_optimal = AUC_conn;

Spatial_amp_optimal = Spatial_amp;
Spatial_conn_optimal = Spatial_conn;


output_eva = output_ortho; 
run('Experiment_evaluate.m')
prob_inclusion_ortho = overall_connectivity;
weighted_inclusion_ortho= overall_connectivity_w;
NRE_conn_ortho = NRE_conn;
NRE_amp_ortho = NRE_amp;
AUC_conn_ortho = AUC_conn;

Spatial_amp_ortho = Spatial_amp;
Spatial_conn_ortho = Spatial_conn;


    save(flnm, 'output_random','output_random_refit','X_random',...
        'location_random','mpp_random','time_record_random',...
        'NRE_conn_random','NRE_amp_random','AUC_conn_random',...
        'output_optimal','output_optimal_refit', 'X_optimal',...
        'location_optimal','mpp_optimal','time_record_optimal',...
        'NRE_conn_optimal','NRE_amp_optimal','AUC_conn_optimal',...
        'output_ortho','output_ortho_refit', 'X_ortho',...
        'location_ortho','mpp_ortho','time_record_ortho',...
        'NRE_conn_ortho','NRE_amp_ortho','AUC_conn_ortho')

end


%%
% covariance of point spread function
A = diag([20, 20, 750]);
for randomseed = 1:10    
    flnm=strcat('../Data/sim-results/A20S', num2str(randomseed),'.mat'); 
    num_sources = 4;  % Number of locations to stimulate in each trial
    run('Parameters_setup_experiment.m')
    % Run analysis and design
    
    design = 0;
    output= struct([]);
    for j = 1:num_threshold
        output(j).alpha = .1*ones(K_z+1,1);
        output(j).mu = zeros(K_z+1,1);
        output(j).s_sq = ones(K_z+1,1);
        %output(j).w_estimate = zeros(K_z+1,1);
    end
    output_refit = output;
    obj_function = @joint_sparsity_weight_entropy; %
    X_g = zeros((N-1)*num_trials_batch + num_trials_first,K_z);
    locations_trials = zeros( (N-1)*num_trials_batch + num_trials_first,num_sources);
    Y_g = zeros((N-1)*num_trials_batch + num_trials_first,num_threshold);
    grid_index = 1:size(pi_dense_all,2);
    
    tic
    locations_record = [];
    time_record = zeros(2*N+1,1);
    time_record(1) = toc;
    
    output_random =cell(N/10,1);
    output_random_refit =cell(N/10,1);
    for outer_i = 1:(N/10)
        t_start = (outer_i-1)*10+1;
        t_end = outer_i*10;
        run('Experiment_full.m')
        output_random{outer_i} = output;
        run('Experiment_refit.m')
        output_random_refit{outer_i} = output_refit;
    end
    X_random = X_g;
    mpp_random = mpp_new;
    time_record_random = time_record;
    location_random = locations_trials;
    
    % Optimal design
    
    design = 1;
    output= struct([]);
    for j = 1:num_threshold
        output(j).alpha = .1*ones(K_z+1,1);
        output(j).mu = zeros(K_z+1,1);
        output(j).s_sq = ones(K_z+1,1);
        %output(j).w_estimate = zeros(K_z+1,1);
    end
    output_refit = output;
    obj_function = @joint_sparsity_weight_entropy; %
    X_g = zeros((N-1)*num_trials_batch + num_trials_first,K_z);
    locations_trials = zeros((N-1)*num_trials_batch + num_trials_first,num_sources);
    Y_g = zeros((N-1)*num_trials_batch + num_trials_first,num_threshold);
    grid_index = 1:size(pi_dense_all,2);
    
    tic
    locations_record = [];
    time_record = zeros(2*N+1,1);
    time_record(1) = toc;
    
    output_optimal =cell(N/10,1);
    output_optimal_refit =cell(N/10,1);
    for outer_i = 1:(N/10)
        t_start = (outer_i-1)*10+1;
        t_end = outer_i*10;
        run('Experiment_full.m')
        output_optimal{outer_i} = output;
        run('Experiment_refit.m')
        output_optimal_refit{outer_i} = output_refit;
    end
    X_optimal = X_g;
    mpp_optimal = mpp_new;
    time_record_optimal = time_record;
    location_optimal = locations_trials;
    
    % near othogonal design
    sigmoid_a=1;
    design = 3;
    output= struct([]);
    for j = 1:num_threshold
        output(j).alpha = .1*ones(K_z+1,1);
        output(j).mu = zeros(K_z+1,1);
        output(j).s_sq = ones(K_z+1,1);
        %output(j).w_estimate = zeros(K_z+1,1);
    end
    output_refit = output;
    obj_function = @joint_sparsity_weight_entropy; %
    X_g = zeros((N-1)*num_trials_batch + num_trials_first,K_z);
    locations_trials = zeros((N-1)*num_trials_batch + num_trials_first,num_sources);
    Y_g = zeros((N-1)*num_trials_batch + num_trials_first,num_threshold);
    grid_index = 1:size(pi_dense_all,2);
    
    tic
    locations_record = [];
    time_record = zeros(2*N+1,1);
    time_record(1) = toc;
    
    output_ortho =cell(N/10,1);
    output_ortho_refit =cell(N/10,1);
    for outer_i = 1:(N/10)
        t_start = (outer_i-1)*10+1;
        t_end = outer_i*10;
        run('Experiment_full.m')
        output_ortho{outer_i} = output;
        run('Experiment_refit.m')
        output_optimal_refit{outer_i} = output_refit;
    end
    X_ortho = X_g;
    mpp_ortho = mpp_new;
    time_record_ortho = time_record;
    location_ortho = locations_trials;
    
    
    % Evaluate the performance: 
% 1, normalized reconstruction error of connectivity
% 2, AUC for connectivity
% 3, normalized amplitudes reconstruction
bandwidth = [200 200];
output_eva = output_random; 
run('Experiment_evaluate.m')

prob_inclusion_random = overall_connectivity;
weighted_inclusion_random = overall_connectivity_w;

NRE_conn_random = NRE_conn;
NRE_amp_random = NRE_amp;
AUC_conn_random = AUC_conn;
Spatial_amp_random = Spatial_amp;
Spatial_conn_random = Spatial_conn;
    

output_eva = output_optimal; 
run('Experiment_evaluate.m')

prob_inclusion_optimal = overall_connectivity;
weighted_inclusion_optimal= overall_connectivity_w;
NRE_conn_optimal = NRE_conn;
NRE_amp_optimal = NRE_amp;
AUC_conn_optimal = AUC_conn;

Spatial_amp_optimal = Spatial_amp;
Spatial_conn_optimal = Spatial_conn;


output_eva = output_ortho; 
run('Experiment_evaluate.m')
prob_inclusion_ortho = overall_connectivity;
weighted_inclusion_ortho= overall_connectivity_w;
NRE_conn_ortho = NRE_conn;
NRE_amp_ortho = NRE_amp;
AUC_conn_ortho = AUC_conn;

Spatial_amp_ortho = Spatial_amp;
Spatial_conn_ortho = Spatial_conn;


    save(flnm, 'output_random','output_random_refit','X_random',...
        'location_random','mpp_random','time_record_random',...
        'NRE_conn_random','NRE_amp_random','AUC_conn_random',...
        'output_optimal','output_optimal_refit', 'X_optimal',...
        'location_optimal','mpp_optimal','time_record_optimal',...
        'NRE_conn_optimal','NRE_amp_optimal','AUC_conn_optimal',...
        'output_ortho','output_ortho_refit', 'X_ortho',...
        'location_ortho','mpp_ortho','time_record_ortho',...
        'NRE_conn_ortho','NRE_amp_ortho','AUC_conn_ortho')

end



%% Computing time:
delta_time_optimal = time_record_optimal(2*(1:N)+1 ) - time_record_optimal(2*(1:N));
delta_time_random = time_record_random(2*(1:N)+1) - time_record_random(2*(1:N));
figure(12)
plot((delta_time_random),'col','r');
hold on;
plot((delta_time_optimal),'col','b');
xlabel('Number of batches');
ylabel('Computing time (s)');
xlim([0,N]);
ylim([-2,4]);
hold off;
%saveas(12,'../Figures/Computing_time.jpg')
%%
figure(1)
histogram([mpp_optimal.amplitudes])
size([mpp_optimal.amplitudes])
figure(2)
histogram([mpp_random.amplitudes])
size([mpp_random.amplitudes])
figure(3)
histogram([mpp_ortho.amplitudes])
size([mpp_ortho.amplitudes])

%%
figure(18)
plot((NRE_conn_random),'col','r');
hold on;
plot((NRE_conn_optimal),'col','b');

%plot((NRE_conn_ortho),'col','g');
xlabel('Number of batches');
ylabel('NRE conn');
xlim([0,N/10+1]);
ylim([0,2]);
hold off;


figure(19)
plot((AUC_conn_random),'col','r');
hold on;
plot((AUC_conn_optimal),'col','b');
%plot((AUC_conn_ortho),'col','g');
xlabel('Number of batches');
ylabel('AUC conn');
xlim([0,N/10+1]);
ylim([0.5,1]);
hold off;

% sum(local_neuron_amplitudes>0)
% TPR v FPR

figure(20)
plot((NRE_amp_random),'col','r');
hold on;
plot((NRE_amp_optimal),'col','b');
%plot((NRE_amp_ortho),'col','g');
xlabel('Number of batches');
ylabel('NRE Amp');
xlim([0,N/10+1]);
ylim([0,2]);
hold off;

figure(21)
plot((Spatial_amp_random),'col','r');
hold on;
plot((Spatial_amp_optimal),'col','b');
%plot((Spatial_amp_ortho),'col','g');
xlabel('Number of batches');
ylabel('Amplitude reconstruction');
xlim([0,N/10+1]);
%ylim([0,2]);
hold off;


figure(22)
plot((Spatial_conn_random),'col','r');
hold on;
plot((Spatial_conn_optimal),'col','b');
%plot((Spatial_conn_ortho),'col','g');
xlabel('Number of batches');
ylabel('Probability reconstruction');
xlim([0,N/10+1]);
%ylim([0,2]);
hold off;



%% Check the design matrix:
XtX_optimal = X_optimal'*X_optimal;
eigs_optimal =eig(XtX_optimal);

XtX_random = X_random'*X_random;
eigs_random =eig(XtX_random);

XtX_ortho = X_ortho'*X_ortho;
eigs_ortho =eig(XtX_ortho);

figure(13)
%histogram(eigs_random)
histogram(diagsq(X_random))

figure(14)
%histogram(eigs_optimal)
histogram(diagsq(X_optimal))


figure(15)
%histogram(eigs_optimal)
histogram(diagsq(X_ortho))

%% Fit a full model to summarize the result
% 
related_mpp_random = mpp_random;

for i = 1:size(X_g,1)
    if size(mpp_random(i).event_times,2) > 0
        indices = mpp_random(i).event_times>evoked_params.stim_start  & mpp_random(i).event_times< (400+evoked_params.stim_start);
        related_mpp_random(i).amplitudes = mpp_random(i).amplitudes(indices);
        related_mpp_random(i).event_times = mpp_random(i).event_times(indices);
    end 
end
amplitude_threshold = quantile([related_mpp_random.amplitudes], (1/num_threshold)*[0:num_threshold]);
amp_counts_random = ones(size(X_g,1),num_threshold);
for j = 1:(num_threshold)
    for i = 1:size(amp_counts_random,1)
        amp_counts_random(i,j) = sum(related_mpp_random(i).amplitudes>amplitude_threshold(j) & related_mpp_random(i).amplitudes<(amplitude_threshold(j+1)+0.01));
    end
end

output_warm = output_random{N/10};
output_post= struct([]);
tic
for j = 1:num_threshold
    fprintf('Running regression on the %d amplitude\n', j);
    Y_n = amp_counts_random(:,j);
    hyperparam_sigma_n =std(Y_n);
    hyperparam_p_connected =  output_warm(j).alpha;
    hyperparam_eta =  output_warm(j).mu;
    hyperparam_sigma_s = ones(size(output_warm(j).mu,1),1);
    options = struct();
    options.alpha =  output_warm(j).alpha;
    options.mu= output_warm(j).mu;
    
    options.verbose= false;
    options.center= 0;
    
    X_temp = X_random(:,:);
    X_temp = [ones(size(X_temp,1),1) X_temp];
    [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs_general(X_temp, Y_n, ...
        hyperparam_sigma_n,hyperparam_sigma_s, hyperparam_p_connected, ...
        hyperparam_eta,options);
    
    output_post(j).alpha = alpha_tmp;
    output_post(j).mu = mu_tmp;
    output_post(j).s_sq = s_sq_tmp;
    output_post(j).w_estimate = alpha_tmp.*mu_tmp;
    
end
toc
output_post_random= output_post;


%%

related_mpp_optimal = mpp_optimal;

for i = 1:size(X_g,1)
    if size(mpp_optimal(i).event_times,2) > 0
        indices = mpp_optimal(i).event_times>evoked_params.stim_start  & mpp_optimal(i).event_times< (400+evoked_params.stim_start);
        related_mpp_optimal(i).amplitudes = mpp_optimal(i).amplitudes(indices);
        related_mpp_optimal(i).event_times = mpp_optimal(i).event_times(indices);
    end 
end
amp_counts_optimal = ones(size(X_g,1),num_threshold);
for j = 1:(num_threshold)
    for i = 1:size(amp_counts_optimal,1)
        amp_counts_optimal(i,j) = sum(related_mpp_optimal(i).amplitudes>amplitude_threshold(j) & related_mpp_optimal(i).amplitudes<(amplitude_threshold(j+1)+0.01));
    end
end

output_warm = output_optimal{N/10};
output_post= struct([]);
tic
for j = 1:(num_threshold)
    fprintf('Running regression on the %d amplitude\n', j);
    Y_n = amp_counts_optimal(:,j);
    hyperparam_sigma_n =std(Y_n);
    hyperparam_p_connected =  output_warm(j).alpha;
    hyperparam_eta =  output_warm(j).mu;
    hyperparam_sigma_s = ones(size(output_warm(j).mu,1),1);
    options = struct();
    options.alpha =  output_warm(j).alpha;
    options.mu= output_warm(j).mu;
    
    options.verbose= false;
    options.center= 0;
    
    X_temp = X_optimal(:,:);
    X_temp = [ones(size(X_temp,1),1) X_temp];
    [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs_general(X_temp, Y_n, ...
        hyperparam_sigma_n,hyperparam_sigma_s, hyperparam_p_connected, ...
        hyperparam_eta,options);
    
    output_post(j).alpha = alpha_tmp;
    output_post(j).mu = mu_tmp;
    output_post(j).s_sq = s_sq_tmp;
    output_post(j).w_estimate = alpha_tmp.*mu_tmp;
    
end
toc
output_post_optimal= output_post;

%% Visualization:
thres = 0.99;
figure(80)

for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        100);
        %neuron_features(i).amplitude(connected_neurons_ind)*35);
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

% for j = 1:(num_threshold)
%     coef = output_random{N/10}(j).alpha(2:end);
%     coef_thres = thres;% quantile(coef,0.98);
%     %sum(coef>coef_thres)
%     
%     potential_neuron_grid = scatter(Z(coef>coef_thres,1), -Z(coef>coef_thres,2),...
%         (amplitude_threshold(j) + amplitude_threshold(j+1)) *35/2, 'filled','o');
%     set(potential_neuron_grid,'MarkerFaceColor','r');
%     alpha(potential_neuron_grid,0.3);
%     hold on
%     
% end
weighted_inclusion_random(weighted_inclusion_random<=0)=0.001;
potential_neuron_grid = scatter(Z(:,1), -Z(:,2), ...
    weighted_inclusion_random*100, 'filled','o');
        
%    (amplitude_threshold(j) + amplitude_threshold(j+1)) *35/2, 'filled','o');
        
    set(potential_neuron_grid,'MarkerFaceColor','r');
    alpha(potential_neuron_grid,0.3);
    hold on
  

temp = scatter(postsyn_position(:,1),...
    -postsyn_position(:,2),...
    182,'filled','o');
set(temp,'MarkerFaceColor','g');
alpha(temp,1);


% Minimize margin and remove ticks
axis off
set(gcf, 'Units','normal')
set(gca, 'Position',[0 0 1 1])
set(gca, 'xtick',[ ])
set(gca, 'ytick',[ ])
view(2)

hold off


figure(90)

for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        100);
        %neuron_features(i).amplitude(connected_neurons_ind)*35);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,0.3);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

%selected_pixels = zeros(num_dense,num_dense,size(amp_related_count_trials,2));
%cent = cell(size(amp_related_count_trials,2),1);

xlim([20,460]);
ylim([-900,-400]);

% for j = 1:(num_threshold)
%     coef = output_optimal{N/10}(j).alpha(2:end);
%     coef_thres = thres;% quantile(coef,0.98);
%     %sum(coef>coef_thres)
%     potential_neuron_grid = scatter(Z(coef>coef_thres,1), -Z(coef>coef_thres,2), ...
%     (amplitude_threshold(j) + amplitude_threshold(j+1)) *35/2, 'filled','o');
%         
%     set(potential_neuron_grid,'MarkerFaceColor','b');
%     alpha(potential_neuron_grid,0.3);
%     hold on
%     
% end

weighted_inclusion_optimal(weighted_inclusion_optimal<=0)=0.001;
potential_neuron_grid = scatter(Z(:,1), -Z(:,2), ...
    weighted_inclusion_optimal*100, 'filled','o');
        
%    (amplitude_threshold(j) + amplitude_threshold(j+1)) *35/2, 'filled','o');
        
    set(potential_neuron_grid,'MarkerFaceColor','b');
    alpha(potential_neuron_grid,0.3);
    hold on
    
temp = scatter(postsyn_position(:,1),...
    -postsyn_position(:,2),...
    182,'filled','o');
set(temp,'MarkerFaceColor','g');
alpha(temp,1);


% Minimize margin and remove ticks
axis off
set(gcf, 'Units','normal')
set(gca, 'Position',[0 0 1 1])
set(gca, 'xtick',[ ])
set(gca, 'ytick',[ ])
view(2)

hold off


%% Visualization:
figure(90)

for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*35);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,0.3);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

xlim([20,460]);
ylim([-900,-400]);

for j = 1:(num_threshold-1)
    coef = output_random_refit{N/10}(j).alpha(2:end);
    coef_thres = 0.99;% quantile(coef,0.98);
    %sum(coef>coef_thres)
    
    potential_neuron_grid = scatter(Z(coef>coef_thres,1), -Z(coef>coef_thres,2),...
        amplitude_threshold(j+1)*35, 'filled','o');
    set(potential_neuron_grid,'MarkerFaceColor','r');
    alpha(potential_neuron_grid,0.3);
    hold on
    
end

for j = 1:(num_threshold-1)
    coef = output_optimal_refit{N/10}(j).alpha(2:end);
    coef_thres = 0.99;% quantile(coef,0.98);
    %sum(coef>coef_thres)
    potential_neuron_grid = scatter(Z(coef>coef_thres,1), -Z(coef>coef_thres,2), ...
        amplitude_threshold(j+1)*35, 'filled','o');
    set(potential_neuron_grid,'MarkerFaceColor','b');
    alpha(potential_neuron_grid,0.3);
    hold on
    
end


temp = scatter(postsyn_position(:,1),...
    -postsyn_position(:,2),...
    182,'filled','o');
set(temp,'MarkerFaceColor','g');
alpha(temp,1);


% Minimize margin and remove ticks
axis off
set(gcf, 'Units','normal')
set(gca, 'Position',[0 0 1 1])
set(gca, 'xtick',[ ])
set(gca, 'ytick',[ ])
view(2)

hold off

%% Draw the stimulated locations:
figure(11)

for i = 1:num_layers
    %connected_neurons_ind = find(neuron_features(i).amplitude);
    connected_neurons_ind = 1:size(neuron_features(i).amplitude,1);
    
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        (neuron_features(i).amplitude(connected_neurons_ind)+1)*35);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,0.3);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

xlim([20,460]);
ylim([-900,-400]);

loc_trials = reshape([location_random( (num_trials_first+1) :end,:) ], [(size(location_random,1)-num_trials_first)*4,1]);
    stimulate = scatter(Z_dense(loc_trials,1),...
        -Z_dense(loc_trials,2),...
        20, 'filled','o');
    set(stimulate,'MarkerFaceColor','r');
   alpha(stimulate,0.1);
    hold off;

    
    
    figure(12)

for i = 1:num_layers
    %connected_neurons_ind = find(neuron_features(i).amplitude);
    connected_neurons_ind = 1:size(neuron_features(i).amplitude,1);
    
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        (neuron_features(i).amplitude(connected_neurons_ind)+1)*35);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,0.3);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

xlim([20,460]);
ylim([-900,-400]);

    loc_trials = reshape([location_optimal( (num_trials_first+1) :end,:) ], [(size(location_optimal,1)-num_trials_first)*4,1]);
    stimulate = scatter(Z_dense(loc_trials,1),...
        -Z_dense(loc_trials,2),...
        20, 'filled','o');
    set(stimulate,'MarkerFaceColor','b');
   alpha(stimulate,0.1);
    hold on;

   
    
    
    figure(13)

for i = 1:num_layers
    %connected_neurons_ind = find(neuron_features(i).amplitude);
    connected_neurons_ind = 1:size(neuron_features(i).amplitude,1);
    
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        (neuron_features(i).amplitude(connected_neurons_ind)+1)*35);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,0.3);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

xlim([20,460]);
ylim([-900,-400]);

    loc_trials = reshape([location_ortho( (num_trials_first+1) :end,:) ], [(size(location_ortho,1)-num_trials_first)*4,1]);
    stimulate = scatter(Z_dense(loc_trials,1),...
        -Z_dense(loc_trials,2),...
        20, 'filled','o');
    set(stimulate,'MarkerFaceColor','g');
   alpha(stimulate,0.1);
    hold on;

%% A map of delta entropy: when picking the first value 

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
       
figure(5)
    entropy_locations(entropy_locations<0) =0; 
    entropies = scatter(Z_dense(:,1),...
        -Z_dense(:,2),...
        abs(entropy_locations)*10+0.01, 'filled','o');
    set(entropies,'MarkerFaceColor','k');
   alpha(entropies,0.4);
    hold on;
    xlim([20,460]);
    ylim([-900,-400]);
    
   [m_id temp_idx] = max(entropy_locations);
    
   peaksfound = scatter(Z_dense(temp_idx,1),...
        -Z_dense(temp_idx,2),...
        100, 'filled','d');
    set(peaksfound,'MarkerFaceColor','r');
   alpha(peaksfound,1);
  
   
   % Minimize margin and remove ticks
   axis off
   set(gcf, 'Units','normal')
   set(gca, 'Position',[0 0 1 1])
   set(gca, 'xtick',[ ])
      set(gca, 'ytick',[ ])
        hold off;
     


 
figure(6)
entropy_locations = entropy_locations - inner_normalized_products(:,temp_idx)*m_id;
    entropy_locations(entropy_locations<0) =0; 
    entropies = scatter(Z_dense(:,1),...
        -Z_dense(:,2),...
        abs(entropy_locations)*10+0.01, 'filled','o');
    set(entropies,'MarkerFaceColor','k');
   alpha(entropies,0.4);
    hold on;
    xlim([20,460]);
    ylim([-900,-400]);
    
   [m_id temp_idx] = max(entropy_locations);
    
   peaksfound = scatter(Z_dense(temp_idx,1),...
        -Z_dense(temp_idx,2),...
        100, 'filled','d');
    set(peaksfound,'MarkerFaceColor','r');
   alpha(peaksfound,1);
  
   
   % Minimize margin and remove ticks
   axis off
   set(gcf, 'Units','normal')
   set(gca, 'Position',[0 0 1 1])
   set(gca, 'xtick',[ ])
      set(gca, 'ytick',[ ])
        hold off;
 
    
        
        
figure(7)
entropy_locations = entropy_locations - inner_normalized_products(:,temp_idx)*m_id;
    entropy_locations(entropy_locations<0) =0; 
    entropies = scatter(Z_dense(:,1),...
        -Z_dense(:,2),...
        abs(entropy_locations)*10+0.01, 'filled','o');
    set(entropies,'MarkerFaceColor','k');
   alpha(entropies,0.4);
    hold on;
    xlim([20,460]);
    ylim([-900,-400]);
    
   [m_id temp_idx] = max(entropy_locations);
    
   peaksfound = scatter(Z_dense(temp_idx,1),...
        -Z_dense(temp_idx,2),...
        100, 'filled','d');
    set(peaksfound,'MarkerFaceColor','r');
   alpha(peaksfound,1);
  
   
   % Minimize margin and remove ticks
   axis off
   set(gcf, 'Units','normal')
   set(gca, 'Position',[0 0 1 1])
   set(gca, 'xtick',[ ])
      set(gca, 'ytick',[ ])
        hold off;
 

saveas(5,'../Figures/Entropy_max_1.jpg')
 
saveas(6,'../Figures/Entropy_max_2.jpg')
      
saveas(7,'../Figures/Entropy_max_3.jpg')
