%% Test the lif-glm model with EM soft assignments


%%
% Loading functions and Data generation
tic
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%%
load('./Environments/l23_cells_lif_glm_fits.mat');
load('./Environments/l23_shapes.mat');

 downsamp=1;
 power_level=l23_cells_first_spike(1).current_data.powers;
 num_power_level=length(power_level);

l23_cells_first_spike(1).stim_start
g=0.01; %membrane time constant [ms]

 current_template=l23_cells_first_spike(1).current_data.current_shape(1:downsamp:end);
 current_template=current_template(100:400); 
 k_offset = 0.002;% l23_cells_first_spike(1).gain; % the spatial mark = k_offset*spatial distance etc
    
%% Paramters in the simulations
num_dense_grid          = [40   40  40  80  80  80  80  80  80  80  80  80  80];
V_homo_grid             = [0    0   0   0   0   0   0   0   0   0   0   0   0];

Aval_grid               = [400  100 900 400 400 400 400 400 400 400 400 400 400]; 
freq_pen_grid           = [1    1   1   1.5 2   1   1   1   1   1   1   1   1];
num_trials_batch_grid   = [20   20  20  20  20  40  100 20  20  20  20  20  20]; 
layer_feature_grid      = [0    0   0   0   0   0   0   1   0   0   0   0   0];
trials_specific_var_grid= [0    0   0   0   0   0   0   0   1   0   0   0   0];
random_prop_grid        = [0    0   0   0   0   0   0   0   0   0.2 0.4 0   0.2];
num_random_grid         = [4    4   4   4   4   4   4   4   4   4   4   2   2]; % number of initial trials

%%
i_setting = 1;
%---------------------------------------------------------------------%
num_dense= num_dense_grid(i_setting); % Number of grids
Aval = Aval_grid(i_setting);
V_homo=V_homo_grid(i_setting);
freq_pen = freq_pen_grid(i_setting);
num_random =num_random_grid(i_setting);
trials_specific_variance= trials_specific_var_grid(i_setting);
random_prop = random_prop_grid(i_setting);
num_trials_batch=num_trials_batch_grid(i_setting);
layer_feature=layer_feature_grid(i_setting);

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
    if V_homo == 1
        cell_feature_priors.Vthre_mean = [0  11  11  11  11  11 11]; % gaussian
        cell_feature_priors.Vreset_mean = [0  -150  -150  -150  -150  -150  -150]; % gaussian
        
    else
        cell_feature_priors.Vthre_mean = [0  9  11  10  13 14 9]; % gaussian
        cell_feature_priors.Vreset_mean = [0  -155  -145  -150  -155  -145  -150]; % gaussian
    end
    cell_feature_priors.Vthre_std = 1* ones(num_layers,1); % gaussian
    cell_feature_priors.Vreset_std =  1* ones(num_layers,1); % gaussian
    %%
    %---------------------------------------------------------------------%
    
    % Parameters for the data generating mechanism
    rng(12242,'twister');
    % load parameters for the model
    
    A = diag([Aval, Aval, 1500]); 
    run('./Data_generation/Parameters_setup_3D.m')
    
    %% Pre-calculation
    % Note: we can use the exp data to estimate pi_dense_all
[pi_dense_all, ~] = get_weights(A,all_locations,all_shape_mat,1e-3,Z_dense);
[pi_dense_local, inner_normalized_products] = get_weights(A,local_locations,local_shape_mat,1e-3,Z_dense);
    %% Visualize  neurons
colrs=cell(2,1);
colrs{1} = [1 0 0];
colrs{2} = [0 0 1];

for i = 1:n_cell_local
h1=scatter3(local_shape_mat{i}(:,1),local_shape_mat{i}(:,2),local_shape_mat{i}(:,3),'.'...
   );
%set(h1, 'MarkerFaceColor',colrs{local_connected(i)+1},'MarkerEdgeColor',colrs{local_connected(i)+1});
%alpha(h1,0.4);

hold on;
end
 xlim([min(Z(:,1))-50, max(Z(:,1))+50]);
    ylim([min(Z(:,2))-50, max(Z(:,2))+50]);
    zlim([min(Z(:,3))-50, max(Z(:,3))+50]);
    
    title('Cells in 3D');
   
    %%
    %num_I_Stim=1; % set the stimulus to be a constant value, e.g., 100
    I_e_vect=current_template;
    evoked_params.stim_start = 1;
    evoked_params.stim_end = length(I_e_vect);
    data_params.T = length(I_e_vect); % total time at 20k Hz
    data_params.dt = 1; % 1/20 ms

    % Stimulation:
    num_sources = 4;  % Number of locations to stimulate in each trial
    grid_index = 1:size(pi_dense_all,2);
    % n_cell = size(pi_dense_all,1);
    % n_cell_local = size(pi_dense_local,1);
    
    k_minimum = 0.001; % minimum stimulus intensity to consider
    
    % Stimulation:
    % Set seed:
    rng(arg1,'twister');
    
    % Parameters
    sqrt_transform = false; % whether to use squared transformation
    % Parameters for the working model
    num_threshold=10; % number of bins to use
    mark = 0; % 0: amplitude; 1: latency.
    obj_function = @joint_sparsity_weight_entropy; %
    
    num_trials_first =max(400, ceil(n_cell_local/num_sources)); % Number of trials in the first batches
    %%
    %---------------------------------------------------------------------%
    % Designing stage
    
    % Design
    %design = 0
    % Run analysis and design
    % Initialize starting values
    output= struct([]);
    for j = 1:num_threshold
        output(j).alpha = .1*ones(n_cell_local*num_power_level+1,1);
        output(j).mu = zeros(n_cell_local*num_power_level+1,1);
        output(j).s_sq = ones(n_cell_local*num_power_level+1,1);
        output(j).threshold = [];
    end
    output_ini = output;
    X_g = zeros(0,n_cell_local*num_power_level);
    locations_trials = zeros(0,num_sources);
    powers_trials= zeros(0,num_sources);
    Y_g = zeros(0,num_threshold);
    % Initialize storage:
    output_crude =cell(N,1);
    time_record = zeros(N,1);
    counts_freq = zeros(size(Z_dense,1)*num_power_level,1);
    % Conducting experiments
    for i_batch = 1:N
       tic
       tstart=toc;
    output_old=output;
    [locations_this_batch, powers_this_batch,counts_freq] = optimal_design_v2(i_batch, num_sources,num_peaks,num_trials_first,num_trials_batch, output, Y_g, ...
      num_power_level,random_prop, counts_freq, pi_dense_local,inner_normalized_products, grid_index, freq_pen, num_samples);
    
    locations_trials = [locations_trials; locations_this_batch];
    powers_trials = [powers_trials; powers_this_batch];

    [mpp_temp, ~, ~] = generate_data(...
        locations_this_batch,powers_this_batch,pi_dense_all,k_offset,k_minimum,all_V_th,all_V_reset,all_gamma,...
        all_amplitudes,all_sigma_across,all_sigma_within, power_level,...
        I_e_vect,g,stoc_mu,stoc_sigma, data_params,bg_params,trials_specific_variance);
    
    [output, Y_g,X_g] = fit_working_model_v2(...
        i_batch,locations_this_batch,powers_this_batch,mpp_temp,Y_g,X_g,output_old,num_threshold,evoked_params,length_memory,k_offset,...
        power_level,mark,pi_dense_local);
    tend=toc;
        
           if i_batch == 1
            mpp= mpp_temp;
        else
            mpp( ((i_batch-2)*num_trials_batch + num_trials_first) + (1:num_trials_batch)) =mpp_temp;
           end
        
        output_crude{i_batch}=output;
        time_record(i_batch)=tend-tstart;
    end
    
%%

crude= cell(N,1);
%threshold =  quantile([mpp.amplitudes], (1/num_threshold)*[0:num_threshold]);
for i_batch = 1:N
    overall_connectivity = zeros(n_cell_local,1);
    overall_mark = zeros(n_cell_local,1);
    for j = 1:num_threshold
        for i_pl = 1:num_power_level
            ind_cell_vec = 1+i_pl+ (0:n_cell_local-1)*num_power_level;
            
        overall_connectivity = max(overall_connectivity, ...
            output_crude{i_batch}(j).alpha(ind_cell_vec));
        for i_cell = 1:n_cell_local
            if overall_connectivity(i_cell) == output_crude{i_batch}(j).alpha(ind_cell_vec(i_cell))
                overall_mark(i_cell)= (output_crude{i_batch}(num_threshold).threshold(j)+...
                    output_crude{i_batch}(num_threshold).threshold(j+1))/2;
            end
        end
        end
    end
    crude{i_batch} = struct;
    crude{i_batch}.mean_gamma = overall_connectivity;
    crude{i_batch}.mean_mu = overall_mark;
    crude{i_batch}.comptime=time_record;
end
%%
jittered_gamma = local_gamma + normrnd(0,0.1,[n_cell_local 1]);
t_seq =[1 2 3 4 5*(1: N/5)];


for ind_t = [1 2 3 10 12]
    
    figure(ind_t)
    %      plot(jittered_gamma,output.EM{ind_t}.mean_gamma,'.','col',[0,1,0,1], 'markersize', 10);
    %      hold on;
    %      plot(jittered_gamma,output.working{ind_t}.mean_gamma,'.','col',[0,0,1,1], 'markersize', 10);
    %      hold on;
    
    plot(jittered_gamma,crude{t_seq(ind_t)}.mean_gamma,'.','col',[1,0,0,1], 'markersize', 10);
    hold on;
    
    line([0 1],[0 1]);
    
    xlim([-0.1,1.1]);
    ylim([-0.1,1.1]);
end

    %%
    %------------------------------------%
    % Estimating the marginal firing rate 
    sigma_unknown=1;
    I_stimuli = I_e_vect;
    T=length(I_e_vect); % total time at 20k Hz
    dt=1;
    t_vect= dt:dt:T;
    if layer_feature == 1
        V_thresholds = local_V_th_layer;
        V_resets = local_V_reset_layer;
    else
        V_thresholds = local_V_th;
        V_resets = local_V_reset;
    end
    
    % Merge the matrix into one 
    
    stimuli_size_local=zeros(size(X_g,1),n_cell_local);
    for i_cell= 1:n_cell_local
        i_cell_vec = (i_cell-1)*num_power_level + (1:num_power_level);
       stimuli_size_local(:,i_cell)= sum(X_g(:,i_cell_vec),2);
    end
    
    n_stimuli_grid=10;
    n_grid_voltage=200;
    t_factor=1;
    
funcs.invlink = @invlink_test;%@(resp) log(1 + exp(resp));%@(x) exp(x);

    [estimated_intensity]=Intensity_v2(stimuli_size_local, n_stimuli_grid,n_grid_voltage,...
        t_vect,t_factor,V_thresholds,V_resets,k_minimum,g,...
        funcs,k_offset,...
        I_stimuli, stoc_mu, stoc_sigma);
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
    %%
    convergence_epsilon = 0.01;
    maxit = 100;
    n_gibbs_sample = 100;
    n_burnin = 100;
    n_skip = 5;
    f_background = bg_params.firing_rate/1000;
    mean_background = bg_params.mean;
    sigma_background  = bg_params.sigma;
    
    
    t_seq = [1 2 3 4 5*(1: N/5)];
    
    output_EM = cell(length(t_seq),1);
    output_Gibbs = cell(length(t_seq),1);
    output_sparse = cell(length(t_seq),1);
    for ind_t = 1:length(t_seq)
        t=t_seq(ind_t); % iterate:
        fprintf('Batch %d\n', t);
        % New data in this batch
        
        idx_trials_this_batch =  1: (num_trials_first + ...
            (t - 1)*num_trials_batch );
        
        n_trial_temp = max(idx_trials_this_batch);
        
       %---------------------------------------------------------------------%
       % Initialization of the full model 
       output= output_crude{t};
       overall_connectivity = zeros(n_cell_local,1);
       overall_mark = zeros(n_cell_local,1);
       for j = 1:num_threshold
           for i_pl = 1:num_power_level
               ind_cell_vec = 1+i_pl+ (0:n_cell_local-1)*num_power_level;
               
               overall_connectivity = max(overall_connectivity, ...
                   output_crude{i_batch}(j).alpha(ind_cell_vec));
               for i_cell = 1:n_cell_local
                   if overall_connectivity(i_cell) == output_crude{i_batch}(j).alpha(ind_cell_vec(i_cell))
                       overall_mark(i_cell)= (output_crude{i_batch}(num_threshold).threshold(j)+...
                           output_crude{i_batch}(num_threshold).threshold(j+1))/2;
                   end
               end
           end
       end
        
        %----------------------------------------------%
        % EM
        sparsity =0;
        [gamma_path mu_path sigma_path total_time]= ...
            EM_fullmodel(mpp(1:n_trial_temp), estimated_intensity(1:n_trial_temp,:),evoked_cell,expected_all, ...
              n_cell_local, overall_connectivity, overall_mark, ones(n_cell_local,1), ...
        convergence_epsilon,f_background, mean_background, sigma_background, sparsity, gamma_threshold,maxit,t_vect);
        %run('./Inference/Simulation_EM_batch.m');
        output_EM{ind_t}.gamma = gamma_path(:,end);
        output_EM{ind_t}.sigma = sigma_path(:,end);
        output_EM{ind_t}.mu = mu_path(:,end);
        output_EM{ind_t}.delta_t = total_time;
        %----------------------------------------------%
        
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