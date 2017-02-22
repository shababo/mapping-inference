%%
clear;
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%%
 flnm=strcat('./Data/example1.mat');
load(flnm);
 flnm=strcat('./Data/truth.mat');
load(flnm);


%% Inference with the working model
num_threshold=10; % number of bins to use
mark = 0; % 0: amplitude; 1: latency.
n_trial = size(stimuli_size,1);
n_cell_local = size(stimuli_size_local,2);

run('./Inference/Simulation_crude.m');

%% Run the full model
%% 1. Estimate the `marginal' firing rate using paramters from LIF-GLM model
% Note: we use the true values in this simulation 
load('./Environments/current_template.mat'); %Contains the vector norm_average_current

I_eg_100mV=abs(norm_average_current(1:20:20*75).*100);
I_eg_50mV=abs(norm_average_current(1:20:20*75).*50);
I_eg_25mV=abs(norm_average_current(1:20:20*75).*25);
I_e=[repmat(I_eg_100mV',1,5) repmat(I_eg_50mV',1,5) repmat(I_eg_25mV',1,5)];
num_I_Stim=1;
I_e_vect=[0;I_e(:,num_I_Stim)];
I_stimuli = I_e_vect;


T=75;
dt=1;
t_vect=0:dt:T;
V_thresholds = local_V_th;
V_resets = local_V_reset;
E_Ls = local_E_L;


        n_stimuli_grid=10;
    
        exact_crossing = 0;
    run('./Inference/Expected_intensity_v2.m');
   %%
    % To improve the speed, we can reduce the computing cost by checking the
   % unique evoked intensity the cell received. 
   
   max(stimuli_size(:,100))
   min(max(stimuli_size))
    % 
    
    % We can further check the similarity of their features:
    

%%
tic
    sigma_unknown=0;

         n_trial_update = 400;
          sigma_unknown=0;
    tstart=toc;
    run('./Inference/Simulation_medium.m');
     tend=toc;
      t_delta = tend-tstart;
       
    flnm=strcat('../../Data/Full_sim/Full_minibatch.mat');
    save(flnm,'t_delta','mpp_new','sigma_samples','gamma_samples','mu_samples','all_connected', 'all_amplitudes',...
        'assignments_samples');
   %%

n_gibbs_sample = 200;
n_burnin = 400;
n_skip = 10;
% Initialize experiment conditions
exact_crossing = 0;
   
   %%
evoked_cell = cell(n_trial,1);
for i_trial = 1:n_trial
    evoked_cell_index = 0; % 0: background evnets
    for i_cell = 1:n_cell_local
        k = stimuli_size_local(i_trial, i_cell);
        if k > 0.01
            evoked_cell_index = [evoked_cell_index i_cell];
        end
    end
    evoked_cell{i_trial} = evoked_cell_index;
end

% 
% soft_assignments_current = cell(n_trial_update,1);
% for i_trial = 1:n_trial_update
%     n_events = length(mpp_new(i_trial).event_times);
%     soft_assignments_current{i_trial,1} = ones(n_events, length(evoked_cell_index))/length(evoked_cell_index);
% end

%%
   tic
         n_trial_update = 800;
          sigma_unknown=0;
    tstart=toc;
    run('./Inference/Simulation_integral.m');
     tend=toc;
      t_delta = tend-tstart;
   
    flnm=strcat('../../Data/Full_sim/Full_minibatch_int.mat');
    save(flnm,'t_delta','sigma_samples','gamma_samples','mu_samples', ...
    'soft_assignments_samples');
   
    %%
    
