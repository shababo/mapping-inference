%%
% Loading functions
tic
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%%  
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
%%
%---------------------------------------------------------------------%

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
    cell_params.shape_gain = local_shape_gain;
    shape_template = l23_cells_for_sim;
[pi_dense_local, inner_normalized_products] = get_weights_v2(cell_params, shape_template,Z_dense);
%%
 cell_params.locations = all_locations;
    cell_params.shape_gain = all_shape_gain;
    shape_template = l23_cells_for_sim;
[pi_dense_all, ~] = get_weights_v2(cell_params, shape_template,Z_dense);
    %% Loading the current template using new template 
    
    load('./Environments/chrome-template-3ms.mat'); 
    
downsamp=1;
power_level=50*ones(5,1);

num_power_level=length(power_level);

 current_template=template(1:downsamp:200); 
  I_e_vect=current_template;
  evoked_params.stim_start = 1;
  evoked_params.stim_end = length(I_e_vect);
 
  
  %num_I_Stim=1; % set the stimulus to be a constant value, e.g., 100
  data_params.T = length(I_e_vect); % total time at 20k Hz
  data_params.dt = 1; % 1/20 ms
%% Setting stimulus parameters:
% Stimulation:
num_sources = 4;  % Number of locations to stimulate in each trial
grid_index = 1:size(pi_dense_all,2);
% n_cell = size(pi_dense_all,1);
% n_cell_local = size(pi_dense_local,1);


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
    % Draw random locations:
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

%%
 
cell_params.V_th = local_V_th;
cell_params.V_reset = local_V_reset;
cell_params.gamma = local_gamma;
cell_params.locations = local_locations;

% The local gains:
cell_params.gain = zeros(n_cell_local,1);
for i_cell = 1 : n_cell_local
    cell_params.gain(i_cell) = l23_cells_for_sim(local_shape_gain(i_cell)).optical_gain;
end

% The local gains:
cell_params.g = zeros(n_cell_local,1);
for i_cell = 1 : n_cell_local
    cell_params.g(i_cell) = l23_cells_for_sim(local_shape_gain(i_cell)).g;
end
%%

    [estimated_intensity]=Intensity_v3(stimuli_size_local, n_stimuli_grid,n_grid_voltage,...
        t_vect,t_factor,k_minimum,...
         cell_params, funcs,...
        I_stimuli, stoc_mu, stoc_sigma);
    %%
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
    emp_all-sum(expected_all,1)'
    %%
    convergence_epsilon = 0.01;
    maxit = 100;
    n_gibbs_sample = 100;
    n_burnin = 100;
    n_skip = 5;
    f_background = bg_params.firing_rate/1000;
    mean_background = bg_params.mean;
    sigma_background  = bg_params.sigma;
    
    
    
    output_EM = cell(1);
    
    %%
%     for ind_t = 1:length(t_seq)
%         t=t_seq(ind_t); % iterate:
%         fprintf('Batch %d\n', t);
%         % New data in this batch
%         
%         idx_trials_this_batch =  1: (num_trials_first + ...
%             (t - 1)*num_trials_batch );
%         
%         n_trial_temp = max(idx_trials_this_batch);
%         
%        %---------------------------------------------------------------------%
%        % Initialization of the full model 
%        output= output_crude{t};
%        overall_connectivity = zeros(n_cell_local,1);
%        overall_mark = zeros(n_cell_local,1);
%        for j = 1:num_threshold
%            for i_pl = 1:num_power_level
%                ind_cell_vec = 1+i_pl+ (0:n_cell_local-1)*num_power_level;
%                
%                overall_connectivity = max(overall_connectivity, ...
%                    output_crude{i_batch}(j).alpha(ind_cell_vec));
%                for i_cell = 1:n_cell_local
%                    if overall_connectivity(i_cell) == output_crude{i_batch}(j).alpha(ind_cell_vec(i_cell))
%                        overall_mark(i_cell)= (output_crude{i_batch}(num_threshold).threshold(j)+...
%                            output_crude{i_batch}(num_threshold).threshold(j+1))/2;
%                    end
%                end
%            end
%        end
        
        %----------------------------------------------%
        % EM
        gamma_old= 0.1*ones(n_cell_local,1);
        mu_old = 2*ones(n_cell_local,1);
        sigma_old = ones(n_cell_local,1);
        sparsity =0;
        [gamma_path mu_path sigma_path total_time]= ...
            EM_fullmodel_v2(mpp(1:n_trial), estimated_intensity(1:n_trial,:),evoked_cell,expected_all, ...
              n_cell_local, gamma_old, mu_old, sigma_old, ...
        convergence_epsilon,f_background, mean_background, sigma_background, sparsity, gamma_threshold,maxit,t_vect);
        %run('./Inference/Simulation_EM_batch.m');
%         output_EM{ind_t}.gamma = gamma_path(:,end);
%         output_EM{ind_t}.sigma = sigma_path(:,end);
%         output_EM{ind_t}.mu = mu_path(:,end);
%         output_EM{ind_t}.delta_t = total_time;
        %----------------------------------------------%
        %%
        gam_est =gamma_path(:,end);
        plot(local_connected+normrnd(0,0.1,[n_cell_local 1] ), gam_est,'.')
         xlim([-0.1,1.1]);
     ylim([-0.1,1.1]);
        %%
        %----------------------------------------------%
        % EM with sparsity  
        sparsity =1;
        [gamma_path mu_path sigma_path total_time]= ...
            EM_fullmodel(mpp(1:n_trial_temp), estimated_intensity(1:n_trial_temp,:),evoked_cell,expected_all, ...
             n_cell_local, overall_connectivity, overall_mark, ones(n_cell_local,1), ...
        convergence_epsilon,f_background, mean_background, sigma_background, sparsity, gamma_threshold,maxit,t_vect);
       output_sparse{ind_t}.gamma = gamma_path(:,end);
        output_sparse{ind_t}.sigma = sigma_path(:,end);
        output_sparse{ind_t}.mu = mu_path(:,end);
        output_sparse{ind_t}.delta_t = total_time;
        %
        
        %--------------------------------------------------------%
        % Gibbs
        gamma_current = overall_connectivity;
        mu_current = overall_mark;
        sigma_current =1*ones(n_cell_local,1);
        %
        mu_m_hyper =overall_mark;
        mu_v_hyper = ones(n_cell_local,1); % inverse of prior variance
        sigma_alpha = ones(n_cell_local,1);
        sigma_beta = ones(n_cell_local,1);
        
        gamma_alpha = ones(n_cell_local,1);
        gamma_beta = ones(n_cell_local,1);
        
          [gamma_path mu_path sigma_path total_time]= ...
            Gibbs_fullmodel(mpp(1:n_trial_temp), estimated_intensity(1:n_trial_temp,:),evoked_cell,expected_all, ...
             n_trial_update, n_cell_local, overall_connectivity, overall_mark, ones(n_cell_local,1), ...
              n_gibbs_sample,n_burnin,n_skip,...
               mu_m_hyper, mu_v_hyper, sigma_alpha, sigma_beta, gamma_alpha, gamma_beta, ...
        convergence_epsilon,f_background, mean_background, sigma_background, sparsity, gamma_threshold,maxit,t_vect);
       output_Gibbs{ind_t}.gamma = gamma_path;
        output_Gibbs{ind_t}.sigma = sigma_path;
        output_Gibbs{ind_t}.mu = mu_path;
        output_Gibbs{ind_t}.delta_t = total_time;
      
    end
    %%
    % record data
    %
    flnm = strcat('./Data/Setting',num2str(i_setting),'Design', num2str(design));
    
    save(strcat(flnm, 'Rep',num2str(arg1),'_Gibbs.mat'),'output_Gibbs');
    save(strcat(flnm, 'Rep',num2str(arg1),'_Sparse.mat'),'output_sparse');
    save(strcat(flnm, 'Rep',num2str(arg1),'_EM.mat'),'output_EM');
    save(strcat(flnm, 'Rep',num2str(arg1),'_Data.mat'),'mpp','stimuli_size_local');
    save(strcat(flnm, 'Rep',num2str(arg1),'_Crude.mat'),'output_crude','time_record','threshold');
    save(strcat(flnm, 'Rep',num2str(arg1),'_Trials.mat'),'locations_trials');
    
    save(strcat(flnm, 'Rep',num2str(arg1),'_Truth.mat'),...
    'bg_params','neuron_features','neuron_locations','num_threshold', ...
        'time_record', 'k_minimum', 'N','num_trials_first','num_trials_batch',...
        'local_locations','local_amplitudes','evoked_params', 'local_V_th','local_V_reset',...
        'local_V_th_layer','local_V_reset_layer',...
        'local_sigma_within','local_sigma_across', 'local_gamma','Z');
    
    
end