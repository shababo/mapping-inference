%%
% Fitting models on the batch data
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%%


flnm = strcat('./Data/Design', num2str(design),'.mat');
load(flnm);
flnm = strcat('./Data/Design', num2str(design),'data.mat');
load(flnm);
flnm = strcat('./Data/Design', num2str(design),'crude.mat');
load(flnm);


%% Some parameters:
n_gibbs_sample = 400;
n_burnin = 800;
n_skip = 20;

lifglm_update = 0;

n_trial_update = 300;

convergence_epsilon = 0.01;
maxit = 100;

sigma_unknown=1;

exact_crossing = 0; % estimate intensities using LIF-GLMs

load('./Environments/current_template.mat'); %Contains the vector norm_average_current
I_eg_100mV=abs(norm_average_current(1:20:20*75).*100);
I_eg_50mV=abs(norm_average_current(1:20:20*75).*50);
I_eg_25mV=abs(norm_average_current(1:20:20*75).*25);
I_e=[repmat(I_eg_100mV',1,5) repmat(I_eg_50mV',1,5) repmat(I_eg_25mV',1,5)];
num_I_Stim=1;
I_e_vect=[0;I_e(:,num_I_Stim)];
I_stimuli = I_e_vect;

% Stochastic components of voltages
stoc_mu=0;stoc_sigma=0.5;
g=0.1; %membrane time constant [ms]

T=75;
dt=1;
t_vect=0:dt:T;

V_thresholds = local_V_th;
V_resets = local_V_reset;

n_stimuli_grid=40;
k_basic = 0.04;

n_stimuli_grid=10;

n_trial = size(stimuli_size_local,1);
n_cell_local = size(stimuli_size_local,2);
%% Prepare the batch data sets

estimated_intensity=cell(0, n_cell_local);
evoked_cell = cell(n_trial,1);

output_EM = cell(N/10,1);
output_Gibbs = cell(N/10,1);

for i_batch_outer = 1:(N/10) % iterate:
    % New data in this batch
    if i_batch_outer == 1
        idx_trials_this_batch = 1:(num_trials_first + num_trials_batch*9);
    else
        idx_trials_this_batch = num_trials_first + ...
            (10*i_batch_outer-10 - 1)*num_trials_batch + 1:(10*num_trials_batch);
    end
    stimuli_temp = stimuli_size_local(idx_trials_this_batch,:);
    run('./Inference/Expected_intensity_v3.m');
    
    estimated_intensity(idx_trials_this_batch,:)=M_intensity;
    for i_trial = idx_trials_this_batch
        evoked_cell_index = 0; % 0: background evnets
        for i_cell = 1:n_cell_local
            k = stimuli_size_local(i_trial, i_cell);
            if k > k_minimum
                evoked_cell_index = [evoked_cell_index i_cell];
            end
        end
        evoked_cell{i_trial} = evoked_cell_index;
    end
    
    n_trial_update = max(idx_trials_this_batch);
    n_trial_temp = max(idx_trials_this_batch);
    
    % Initialize the estimates
    output = output_random{i_batch_outer};
    
    overall_connectivity = zeros(n_cell_local,1);
    overall_connectivity_w = zeros(n_cell_local,1);
    
    overall_mark = zeros(n_cell_local,1);
    normalized_constants = zeros(n_cell_local,1);
    for j = 1:num_threshold
        overall_connectivity = max(overall_connectivity, ...
            output(j).alpha(2:end));
        for i_cell = 1:n_cell_local
            if overall_connectivity(i_cell) == output(j).alpha(i_cell+1)
                overall_mark(i_cell)= (threshold(j)+threshold(j+1))/2;
            end
        end
        
    end
    
    gamma_old = overall_connectivity;
    mu_old = overall_mark;
    sigma_old = 1*ones(n_cell_local,1);
    f_background = bg_params.firing_rate/1000;
    w_background = 1/range([mpp_new(1:n_trial_temp).amplitudes]);
    
    expected_all = zeros(n_trial_temp,i_cell);
    for i_trial = 1:n_trial_temp
        for i_cell = 1:n_cell_local
            expected_all(i_trial,i_cell)=sum(estimated_intensity{i_trial, i_cell} );
        end
    end
    
    run('./Inference/Simulation_EM_batch.m');
    output_EM{i_batch_outer}.sigma_samples = sigma_samples;
    output_EM{i_batch_outer}.gamma_samples = gamma_samples;
    output_EM{i_batch_outer}.mu_samples = mu_samples;
    output_EM{i_batch_outer}.soft_assignments_samples= soft_assignments_samples;
    
    
    gamma_current = overall_connectivity;
    mu_current = overall_mark;
    sigma_current =4*ones(n_cell_local,1);
    %
    mu_m_hyper =overall_mark;
    mu_v_hyper = ones(n_cell_local,1); % inverse of prior variance
    sigma_alpha = ones(n_cell_local,1);
    sigma_beta = ones(n_cell_local,1);
    
    gamma_alpha = ones(n_cell_local,1);
    gamma_beta = ones(n_cell_local,1);
    
    f_background = bg_params.firing_rate/1000;
    w_background = 1/range([mpp_new(1:n_trial_temp).amplitudes]); % size distribution for background events
    
    
    run('./Inference/Simulation_integral_batch.m');
    output_Gibbs{i_batch_outer}.sigma_samples = sigma_samples;
    output_Gibbs{i_batch_outer}.gamma_samples = gamma_samples;
    output_Gibbs{i_batch_outer}.mu_samples = mu_samples;
    output_Gibbs{i_batch_outer}.soft_assignments_samples= soft_assignments_samples;
    
end


flnm=strcat('./Data/Design',num2str(design), 'Batch_EM.mat');
save(flnm,'output_EM');
flnm=strcat('./Data/Design',num2str(design),'Batch_Gibbs.mat');
save(flnm,'output_Gibbs');
flnm=strcat('./Data/Design',num2str(design),'Batch_Intensity.mat');
save(flnm,'estimated_intensity');

