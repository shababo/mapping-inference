%%
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Paramters in the simulations
num_dense_grid          = [40];
freq_pen_grid           = [1];
num_trials_batch_grid   = [20];
trials_specific_var_grid= [0];
random_prop_grid        = [1];
num_random_grid         = [4];
i_setting = 1;

num_dense= num_dense_grid(i_setting); % Number of grids
freq_pen = freq_pen_grid(i_setting);
num_trials_batch = num_trials_batch_grid(i_setting);
trials_specific_variance= trials_specific_var_grid(i_setting);
random_prop = random_prop_grid(i_setting);
num_trials_batch=num_trials_batch_grid(i_setting);
%% Loading templates from real data
load('./Environments/l23_cells_for_sim.mat');

num_types_cell = length(l23_cells_for_sim);
% normalized the cell shapes
for i = 1:num_types_cell
    temp=l23_cells_for_sim(i).shape;
    temp_max = max(max(max(temp)));
    l23_cells_for_sim(i).shape = temp/temp_max;
end
% Manually set the 0 gain to 0.01 (to avoid connected cell that are not
% activated)
l23_cells_for_sim(13).optical_gain=0.01;
%

%% Load templates (inference)


load('./Environments/l23_template_cell.mat');
%l23_average_shape
temp=l23_average_shape;
temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;

%%
n_trial_update = 200;

length_memory = num_trials_batch;
num_peaks = 2*num_trials_batch;
gamma_threshold = 0.1; % For the sparse EM algorithm

num_samples=50; % Number of samples to estimate the expected entropies
num_sources = 4;  % Number of locations to stimulate in each trial
N=4000/num_trials_batch; % Number of batches

cell_feature_priors=struct();
num_layers = 7;
cell_feature_priors.Vthre_mean = [0  15  15  15  15  15 15]; % gaussian
cell_feature_priors.Vreset_mean = [0  -1e3  -1e3  -1e3  -1e3  -1e3  -1e3]; % gaussian
cell_feature_priors.Vthre_std = 0* ones(num_layers,1); % same
cell_feature_priors.Vreset_std =  0* ones(num_layers,1);
%%
%---------------------------------------------------------------------%
% Parameters for the data generating mechanism
rng(12242,'twister');
% load parameters for the model
run('./Data_generation/Parameters_setup_3D.m')
%% Pre-calculation
% Note: use the standard template for inference when testing robustness
cell_params.locations = local_locations;
cell_params.shape_gain = ones(n_cell_local,1);
shape_template = struct();
shape_template.shape= l23_average_shape;

[pi_dense_local, inner_normalized_products] = get_weights_v2(cell_params, shape_template,Z_dense);
%%
cell_params.locations = all_locations;
cell_params.shape_gain = all_shape_gain;
shape_template = l23_cells_for_sim;
[pi_dense_all, ~] = get_weights_v2(cell_params, shape_template,Z_dense);
%% Loading the current template using new template
load('./Environments/chrome-template-3ms.mat');
downsamp=1;
power_level=[20 30 40 50 60];% arbitrary numbers to make it work
num_power_level=length(power_level);
current_template=template(1:downsamp:200);
I_e_vect=current_template;
evoked_params.stim_start = 1;
evoked_params.stim_end = length(I_e_vect);
data_params.T = length(I_e_vect); % total time at 20k Hz
data_params.dt = 1; % 1/20 ms
%% Setting stimulus parameters:
% Stimulation:
num_sources = 4;  % Number of locations to stimulate in each trial
grid_index = 1:size(pi_dense_all,2);

% Parameters
sqrt_transform = false; % whether to use squared transformation
% Parameters for the working model
num_threshold=10; % number of bins to use
mark = 0; % 0: amplitude; 1: latency.
obj_function = @joint_sparsity_weight_entropy; %
num_trials_first =max(400, ceil(n_cell_local/num_sources)); % Number of trials in the first batches
k_minimum = 0.001; % minimum stimulus intensity to consider
%%
%---------------------------------------------------------------------%
% Design stage
% initialization
output= struct([]);
for j = 1:num_threshold
    output(j).alpha = .1*ones(n_cell_local*num_power_level+1,1);
    output(j).mu = zeros(n_cell_local*num_power_level+1,1);
    output(j).s_sq = ones(n_cell_local*num_power_level+1,1);
    output(j).threshold = [];
end
X_g = zeros(0,n_cell_local*num_power_level);
locations_trials = zeros(0,num_sources);
powers_trials= zeros(0,num_sources);
Y_g = zeros(0,num_threshold);
counts_freq = zeros(size(Z_dense,1)*num_power_level,1);
%%
cell_params.V_th = all_V_th;
cell_params.V_reset = all_V_reset;
cell_params.gamma = all_gamma;
cell_params.amplitudes = all_amplitudes;
cell_params.sigma_across = all_sigma_across;
cell_params.sigma_within = all_sigma_within;
cell_params.locations = all_locations;
cell_params.shape_gain = all_shape_gain;

stoc_params.mu=stoc_mu;
stoc_params.sigma=stoc_sigma;
%%
for i_batch= 1:50
    tic
    tstart=toc;
    output_old=output;
    [locations_this_batch, powers_this_batch,counts_freq] = optimal_design_v2(i_batch, num_sources,num_peaks,num_trials_first,num_trials_batch, output, Y_g, ...
        num_power_level,random_prop, counts_freq, pi_dense_local,inner_normalized_products, grid_index, freq_pen, num_samples);
    
    locations_trials = [locations_trials; locations_this_batch];
    powers_trials = [powers_trials; powers_this_batch];
    
    [mpp_temp, ~, ~] = generate_data_v2(...
        locations_this_batch,powers_this_batch,pi_dense_all,k_minimum,cell_params, shape_template, power_level,...
        I_e_vect,stoc_params, data_params,bg_params,trials_specific_variance);
    
    if i_batch == 1
        mpp= mpp_temp;
    else
        mpp( ((i_batch-2)*num_trials_batch + num_trials_first) + (1:num_trials_batch)) =mpp_temp;
    end
end
%%
%------------------------------------%
% Estimating the marginal firing rate
sigma_unknown=1;
I_stimuli = I_e_vect;
T=length(I_e_vect); % total time at 20k Hz
dt=1;
t_vect= dt:dt:T;

stimuli_size_local=zeros(length(mpp),n_cell_local);
for l = 1:length(mpp)
    for m = 1:size(locations_trials,2)
        stimuli_size_local(l,:)  = stimuli_size_local(l,:)+( pi_dense_local(:,locations_trials(l,m)).*power_level(powers_trials(l,m)))';
    end
end

n_stimuli_grid=10;
n_grid_voltage=400;
t_factor=1;


funcs.invlink = @invlink_test;%@(resp) log(1 + exp(resp));%@(x) exp(x);

n_trial = size(stimuli_size_local,1);
evoked_cell = cell(n_trial,1);
for i_trial = 1:n_trial
    evoked_cell_index = 0; % 0: background evnets
    for i_cell = 1:n_cell_local
        k = stimuli_size_local(i_trial, i_cell);
        if k > k_minimum
            evoked_cell_index = [evoked_cell_index i_cell];
        end
    end
    evoked_cell{i_trial} = evoked_cell_index;
end


%%
cell_params.V_th = local_V_th;
cell_params.V_reset = local_V_reset;
cell_params.locations = local_locations;

% The local gains:
cell_params.gain = zeros(n_cell_local,1);
for i_cell = 1 : n_cell_local
    %cell_params.gain(i_cell) = l23_cells_for_sim(local_shape_gain(i_cell)).optical_gain;
    cell_params.gain(i_cell) = mean([l23_cells_for_sim.optical_gain]);
end

% The local g:
cell_params.g = zeros(n_cell_local,1);
for i_cell = 1 : n_cell_local
    %cell_params.g(i_cell) = l23_cells_for_sim(local_shape_gain(i_cell)).g;
    cell_params.g(i_cell) =  mean([l23_cells_for_sim.g]);
end
%%
[estimated_intensity]=Intensity_v4(stimuli_size_local, n_stimuli_grid,n_grid_voltage,...
    t_vect,t_factor,k_minimum,...
    cell_params, funcs,...
    I_stimuli, stoc_mu, stoc_sigma);

expected_all = zeros(n_trial,n_cell_local);
for i_trial = 1:n_trial
    for i_cell = 1:n_cell_local
        expected_all(i_trial,i_cell)=sum(estimated_intensity{i_trial, i_cell} );
    end
end
%---------------------------------------------------------------------%
%% Sanity check: expected counts v.s. empirical counts
emp_all = zeros(n_cell_local,1);
for i_cell = 1:n_cell_local
    emp_all(i_cell) = sum([mpp.assignments] ==local_index(i_cell));
end
emp_all(local_connected)-sum(expected_all(:,local_connected),1)'
%%
convergence_epsilon = 0.01;
maxit = 100;
n_gibbs_sample = 100;
n_burnin = 100;
n_skip = 5;
f_background = bg_params.firing_rate/1000;
mean_background = bg_params.mean;
sigma_background  = bg_params.sigma;

gain_old = 0.003*ones(n_cell_local,1);
gamma_old= 0.009*ones(n_cell_local,1);
mu_old = 2*ones(n_cell_local,1);
sigma_old = ones(n_cell_local,1);
sparsity =0;

[gamma_path mu_path sigma_path total_time soft_assignments]= ...
    EM_fullmodel_v2(mpp(1:n_trial), estimated_intensity(1:n_trial,:),evoked_cell,expected_all, ...
    n_cell_local, gamma_old, mu_old, sigma_old, ...
    convergence_epsilon,f_background, mean_background, sigma_background, sparsity, gamma_threshold,maxit,t_vect);

gamma_current= gamma_path(:,end);
sigma_current = sigma_path(:,end);
mu_current = mu_path(:,end);
%gain_current = median(gains_sample,2);
%%
gamma_current
mpp.assignments
%% Remove the disconnected cells 
connected_threshold = 1e-4;
connected_cells = gamma_current > connected_threshold;

%% 
stimuli_size_temp = stimuli_size_local(:,connected_cells);
n_cell_temp = sum(connected_cells);

cell_params.V_th = local_V_th(connected_cells);
cell_params.V_reset = local_V_reset(connected_cells);
cell_params.locations = local_locations;
% The local gains:
cell_params.gain = zeros(n_cell_temp,1);
for i_cell = 1 : n_cell_temp
    %cell_params.gain(i_cell) = l23_cells_for_sim(local_shape_gain(i_cell)).optical_gain;
    cell_params.gain(i_cell) = mean([l23_cells_for_sim.optical_gain]);
end

% The local g:
cell_params.g = zeros(n_cell_temp,1);
for i_cell = 1 : n_cell_temp
    %cell_params.g(i_cell) = l23_cells_for_sim(local_shape_gain(i_cell)).g;
    cell_params.g(i_cell) =  mean([l23_cells_for_sim.g]);
end
%
[estimated_intensity]=Intensity_v4(stimuli_size_temp, n_stimuli_grid,n_grid_voltage,...
    t_vect,t_factor,k_minimum,...
    cell_params, funcs,...
    I_stimuli, stoc_mu, stoc_sigma);

expected_all = zeros(n_trial,n_cell_temp);
for i_trial = 1:n_trial
    for i_cell = 1:n_cell_temp
        expected_all(i_trial,i_cell)=sum(estimated_intensity{i_trial, i_cell} );
    end
end


evoked_cell_temp = cell(n_trial,1);
for i_trial = 1:n_trial
    evoked_cell_index = 0; % 0: background evnets
    for i_cell = 1:n_cell_temp
        k = stimuli_size_temp(i_trial, i_cell);
        if k > k_minimum
            evoked_cell_index = [evoked_cell_index i_cell];
        end
    end
    evoked_cell_temp{i_trial} = evoked_cell_index;
end



gain_old = 0.003*ones(n_cell_temp,1);
gamma_old= 0.009*ones(n_cell_temp,1);
mu_old = 2*ones(n_cell_temp,1);
sigma_old = ones(n_cell_temp,1);

%% Start of the iteration:
normalized_change_outer = convergence_epsilon_outer + 1;
num_iter = 0;
convergence_epsilon_outer= 0.001;

while (normalized_change_outer > convergence_epsilon_outer) & (num_iter < maxit)
    num_iter = num_iter+1;
    %----------------------------------------------%
    % EM
    [gamma_path mu_path sigma_path total_time soft_assignments]= ...
        EM_fullmodel_v2(mpp(1:n_trial), estimated_intensity(1:n_trial,:),evoked_cell_temp,expected_all, ...
        n_cell_local, gamma_old, mu_old, sigma_old, ...
        convergence_epsilon,f_background, mean_background, sigma_background, sparsity, ...
        gamma_threshold,maxit,t_vect);
    
    % lif-glm updates:
    % Use Monte-Carlo method to update lif-glm parameters based on soft
    % assignments
    % should be turned into a function
    num_MC_lifglm = 20;
    
    gains_sample = zeros(n_cell_temp,num_MC_lifglm);
    for i_MC = 1:num_MC_lifglm
        
        % Draw hard assignments using the soft ones:
        mpp_rand = mpp;
        trials_evoked_cell = cell(n_cell_temp,1);
        for i_trial = 1:n_trial
            n_event = length(mpp_rand(i_trial).event_times);
            cell_list = evoked_cell_temp{i_trial};
            mpp_rand(i_trial).assignments(:)=0;
            for i_event = 1:n_event
                this_assignments = soft_assignments{i_trial}(i_event,:);
                r_temp = rand(1);
                i_cell = min(find(r_temp<cumsum(this_assignments)));
                mpp_rand(i_trial).assignments(i_event) = cell_list(i_cell);
                if cell_list(i_cell) ~=0
                    trials_evoked_cell{cell_list(i_cell)} =  [trials_evoked_cell{cell_list(i_cell)} i_trial];
                end
            end
        end
        
        % Reformat the hard threshold for fitting lif glm
        cell_data = cell(n_cell_temp,1);
        for i_cell = 1:n_cell_temp
            trial_list = trials_evoked_cell{i_cell};
            N_cell = length(trial_list);
            if N_cell>0
                cell_data{i_cell}=struct();
                cell_data{i_cell}.responses = zeros(N_cell, length(I_e_vect));
                cell_data{i_cell}.stims = zeros(N_cell, length(I_e_vect));
                for idx_trial = 1:N_cell
                    i_trial = trial_list(idx_trial);
                    event_times = round(mpp_rand(i_trial).event_times(mpp_rand(i_trial).assignments == i_cell));
                    cell_data{i_cell}.responses(idx_trial,event_times) = 1;
                    cell_data{i_cell}.stims(idx_trial,:) =  I_e_vect*stimuli_size_temp(i_trial, i_cell); 
                end
            end
        end
         lif_glm_gains= zeros(n_cell_temp,1);
        for i_cell = 1:n_cell_temp
            trial_list = trials_evoked_cell{i_cell};
            N_cell = length(trial_list);
            in_params.g = cell_params.g(i_cell);
            if N_cell>2
                stims_ind = ones(size(cell_data{i_cell}.stims,1),1);
                % LIF-GLM fits
                %-------------------------------------%
                [stats_conv] = fit_lifglm_v2( cell_data{i_cell}.responses, cell_data{i_cell}.stims,stims_ind,in_params);
                %-------------------------------------%
                lif_glm_gains(i_cell)=stats_conv.beta(2);
            else
                lif_glm_gains(i_cell)=cell_params.gain(i_cell);
            end
        end
        
        gains_sample(:,i_MC)=lif_glm_gains;
    end
    
    % Evaluate the updates:
    gamma_current= gamma_path(:,end);
    sigma_current = sigma_path(:,end);
    mu_current = mu_path(:,end);
    gain_current = median(gains_sample,2);
    
    normalized_change_outer = norm(gamma_current - gamma_old)/(norm(gamma_old)+1) + norm(mu_current - mu_old)/(norm(mu_old)+1)+...
        norm(sigma_current - sigma_old)/(norm(sigma_old)+1)+norm(gain_current-gain_old)/(norm(gain_old)+1);
    
    gamma_old =  gamma_current;
    sigma_old = sigma_current;
    mu_old = mu_current;
    gain_old = gain_current;
    
    fprintf('Changes %d\n', normalized_change_outer);
    % Update the intensities
cell_params.V_th = local_V_th(connected_cells);
cell_params.V_reset = local_V_reset(connected_cells);
cell_params.locations = local_locations;
    % The local gains: update gains based on the lif-glm fits
    cell_params.gain = gain_current;
    [estimated_intensity]=Intensity_v4(stimuli_size_temp, n_stimuli_grid,n_grid_voltage,...
        t_vect,t_factor,k_minimum,...
        cell_params, funcs,...
        I_stimuli, stoc_mu, stoc_sigma);
    
    expected_all = zeros(n_trial,n_cell_temp);
    for i_trial = 1:n_trial
        for i_cell = 1:n_cell_temp
            expected_all(i_trial,i_cell)=sum(estimated_intensity{i_trial, i_cell} );
        end
    end     
end
%%
gam_est =gamma_path(:,end);
plot(local_gamma+normrnd(0,0.1,[n_cell_local 1] ), gam_est,'.','markers',12)
xlim([-0.1,1.1]);
ylim([-0.1,1.1]);
%%

plot(sqrt(cell_params.gain(connected_cells))+normrnd(0,0.01,[n_cell_temp 1] ), gain_current,'.','markers',12)
xlim([-0.05,0.2]);
ylim([-0.05,0.2]);