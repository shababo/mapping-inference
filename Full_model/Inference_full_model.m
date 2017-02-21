%%
clear;
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%%
 flnm=strcat('./Data/example1.mat');
load(flnm);


%% Inference with the working model
num_threshold=10; % number of bins to use
mark = 0; % 0: amplitude; 1: latency.
run('./Inference/Simulation_crude.m');
   

%% Run the full model
%% 1. Estimate the `marginal' firing rate using paramters from LIF-GLM model
% Note: we use the true values in this simulation 

I_stimuli = I_e_vect;
V_thresholds = all_V_th;
V_resets = all_V_reset;
E_Ls = all_E_L;

        n_stimuli_grid=10;
    
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
    save(flnm,'t_delta','mpp_new','sigma_samples','gamma_samples','mu_samples','all_connected', 'all_amplitudes', ...
    'soft_assignments_samples');
   
    %%
    
