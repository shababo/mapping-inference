%% Loading functions and Data generation
clear;
addpath(genpath('../psc-detection'),genpath('../mapping-inference'),genpath('../mapping-core'));
%% Parameters for the data generating mechanism
rng(12242,'twister');

% Parameters on the grid and experimental design 
region_width = 500;
region_height = 500;

% load parameters for the model
run('Parameters_setup_ground_truth.m')

%% 
% Pre-processing for  experiment
num_samples=50; % Number of samples to estimate the expected entropies
num_trials_batch=20;

num_sources = 4;  % Number of locations to stimulate in each trial
num_peaks = 30;
% We consider new locations on a grid (rather than the cell somas) in order
% to gain extra information (e.g., to be able to distinguish two
% neighbouring cells)

num_dense=40; % Grid density


run('Parameters_setup_experiment.m')


%% Run analysis and design 
N=200; % Number of batches
N_r = 50; % Number of batches with random designs 

num_threshold=20;
amp_max = 15;
amp_min = 0;
amplitude_threshold = (0:num_threshold) * (amp_max-amp_min)/num_threshold;

run('Experiment_full.m')
%% Computing time:
delta_time = time_record(2: (N)) - time_record(1:(N-1));
figure(12)
plot(delta_time);
hold on;
xlabel('Number of batches');
ylabel('Computing time (s)');
xlim([0,N]);
ylim([0,10]);
hold off;

%% Fit a semi-full model to summarize the result
output_warm = output;
output_post= struct([]);
tic
for j = 1:size(amp_related_count_trials,2)
    Y_n = Y_g(:,j);
    hyperparam_sigma_n =std(Y_n);
    hyperparam_p_connected =  output_warm(j).alpha;
    hyperparam_eta =  output_warm(j).mu;
    hyperparam_sigma_s = ones(size(output_warm(j).mu,1),1);
    options = struct();
    options.alpha =  output_warm(j).alpha;
    options.mu= output_warm(j).mu;
    
    options.verbose= false;
    options.center= 0;
    
    X_temp = X_g(:,:);
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

for j = 1:(num_threshold-1)
    coef = output(j).alpha(2:end);
    coef_thres = 0.8;% quantile(coef,0.98);
    potential_neuron_grid = scatter(Z(coef>coef_thres,1), -Z(coef>coef_thres,2), amplitude_threshold(j+1)*25, 'filled','o');
    set(potential_neuron_grid,'MarkerFaceColor','r');
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