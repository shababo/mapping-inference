%%
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Loading templates from real data
load('./Environments/l23_cells_for_sim.mat');
% We only use the main gain in this analysis 
mean_gain = mean([l23_cells_for_sim.optical_gain])/2;
%% Load templates for inference 
load('./Environments/l23_template_cell.mat');
%l23_average_shape
temp=l23_average_shape;
temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;
%% Set seed for reproducibility 
rng(12242,'twister');
% load real data
%load('./Environments/05082017_s3c1_t345_mpp_stim_data.mat') %for nuc_locs and stims 
load('./Environments/5_27_s1c1_mpp_and_stim_data.mat') 

% %% Draw some plots..
% 
% figure(1)
% temp2 = scatter(cell_locs(:,2),cell_locs(:,1),...
%     50);
% set(temp2,'MarkerFaceColor','r');
% alpha(temp2,0.5);
% 
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])
% 
% axis off;
% 
% %%
% figure(2)
% temp2 = scatter(target_locs(:,2),target_locs(:,1),...
%     50);
% set(temp2,'MarkerFaceColor','r');
% alpha(temp2,0.5);
% 
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])
% 
% axis off;

%% Pre-calculation
% Note: use the standard template for inference when testing robustness
cell_params.locations = cell_locs;
local_locations = cell_locs;


n_cell_local = size(cell_locs,1);
cell_params.shape_gain = ones(n_cell_local,1);
shape_template = struct();
shape_template.shape= l23_average_shape;
Z_dense =target_locs;
[pi_dense_local, inner_normalized_products] = get_weights_v2(cell_params, shape_template,Z_dense);
%% Loading the current template using new template
load('./Environments/chrome-template-3ms.mat');
max_time = 300;
downsamp=1;
current_template=template(1:downsamp:300);
I_e_vect=current_template;
evoked_params.stim_start = 1;
evoked_params.stim_end = length(I_e_vect);
data_params.T = length(I_e_vect); % total time at 20k Hz
data_params.dt = 1; % 1/20 ms
k_minimum = 0.001; % minimum stimulus intensity to consider

%%
%------------------------------------%
% Estimating the marginal firing rate
sigma_unknown=1;
I_stimuli = I_e_vect;
T=length(I_e_vect); % total time at 20k Hz
dt=1;t_vect= dt:dt:T;

sd_range=1.5;
powers_trials = stim_pow;
n_trial = size(stim_pow,1);

stimuli_size_local=zeros(n_trial,n_cell_local);
% locations_trials = arrayfun(@(x) find(arrayfun(@(y) isequal(y,x),nuc_locs)),Z_dense);
locations_trials=target_inds;

%locations_trials  = locations_trials';
for l = 1:n_trial
    for m = 1:size(locations_trials,2)
        if isnan(locations_trials(l,m))
            
        else
            stimuli_size_local(l,:)  = stimuli_size_local(l,:)+( pi_dense_local(:,locations_trials(l,m)).*powers_trials(l))';
        end
       
    end
end

n_stimuli_grid=40;
n_grid_voltage=400;
t_factor=1;
funcs.invlink = @invlink_test;%@(resp) log(1 + exp(resp));%@(x) exp(x);

%% Removing cells that are rarely stimulated
stim_threshold = 20;
stimulated_cells = sum(stimuli_size_local>stim_threshold )>10;
% also throw away cells with Z>140
stimulated_cells = stimulated_cells & (local_locations(:,3)<140)';

sum(stimulated_cells)

%%

stimuli_size_temp = stimuli_size_local(:,stimulated_cells);
n_cell_temp = sum(stimulated_cells);

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

v_reset_known = -4e3;
v_th_known=15;
% The local gains:
cell_params.gain = zeros(n_cell_local,1);
cell_params.g = zeros(n_cell_local,1);

% delay_params.mean=1.75*20;
% delay_params.std=2.7*20;

delay_params.mean=0;
delay_params.std=0;

cell_params.V_th =  v_th_known*ones(sum(stimulated_cells),1);
cell_params.V_reset = -v_reset_known*ones(sum(stimulated_cells),1);

cell_params.gain = zeros(n_cell_temp,1);
cell_params.g = zeros(n_cell_temp,1);

for i_cell = 1 : n_cell_temp
    %cell_params.gain(i_cell) = l23_cells_for_sim(local_shape_gain(i_cell)).optical_gain;
    cell_params.gain(i_cell) = mean_gain;
    cell_params.gain_sd(i_cell)= std([l23_cells_for_sim.optical_gain]);
    %cell_params.g(i_cell) = l23_cells_for_sim(local_shape_gain(i_cell)).g;
    cell_params.g(i_cell) =  mean([l23_cells_for_sim.g]);
end
%

%%
min_time = 60;
mpp_copy=mpp;
for i_trial = 1:n_trial
    if mpp(i_trial).num_events >0
        range_idx = mpp(i_trial).times<max_time & mpp(i_trial).times>min_time;
        
        mpp(i_trial).num_events = sum(range_idx);
        mpp(i_trial).times = mpp(i_trial).times(range_idx);
        mpp(i_trial).amp = mpp(i_trial).amp(range_idx);
    end
    
end

%%
trial_counts = 0;
event_counts = 0;
for i_trial = 1:n_trial
    %if stim_pow(i_trial)==50
       trial_counts = trial_counts+1;
       event_counts = event_counts +sum(mpp(i_trial).times>60 & mpp(i_trial).times<100);
    %end
end
background_rate = event_counts/trial_counts/40;
mpp_copy=mpp;
%%
[event_rates,expected_all] =Intensity_v6(stimuli_size_temp,  n_stimuli_grid,n_grid_voltage,...
    t_vect,t_factor,k_minimum,...
    cell_params, funcs,...
    I_stimuli,sd_range,delay_params,...
    mpp);

%---------------------------------------------------------------------%
%%

background_update=0;

f_background = background_rate;


mean_background = 0; %not actually used
sigma_background  =1;%not actually used



gain_old = 0.003*ones(n_cell_temp,1);
gamma_old= 0.009*ones(n_cell_temp,1);

mu_old = 2*ones(n_cell_temp,1);
sigma_old = ones(n_cell_temp,1);

sparsity =0; %use sparsity or not 
gamma_threshold = 0.1;
use_size =0;
convergence_epsilon = 0.01;
maxit = 100;


[gamma_path mu_path sigma_path total_time soft_assignments bg_rate]= ...
    EM_fullmodel_v3(mpp(1:n_trial), ...
    event_rates,...
    evoked_cell_temp,expected_all, ...
    n_cell_local, gamma_old, mu_old, sigma_old, ...
    convergence_epsilon,f_background, mean_background, sigma_background, ...
    sparsity, gamma_threshold,maxit,t_vect,use_size,background_update);

gamma_ini= gamma_path(:,end);
sigma_ini = sigma_path(:,end);
mu_ini = mu_path(:,end);
%gain_current = median(gains_sample,2);

%% Select cells:
% connected_cells = gamma_ini>0; %keeping all the cells 
% 
% n_cell_temp=sum(connected_cells);
% cell_params.locations = cell_locs(connected_cells,:);
% cell_params.shape_gain = ones(n_cell_temp,1);
% shape_template = struct();
% shape_template.shape= l23_average_shape;
% [pi_dense_local, inner_normalized_products] = get_weights_v2(cell_params, shape_template,Z_dense);
% 
%  %%
% 
% %locations_trials  = locations_trials';
% stimuli_size_local= zeros(n_trial,sum(connected_cells));
% for l = 1:n_trial
%     for m = 1:size(locations_trials,2)
%         if isnan(locations_trials(l,m))
%             
%         else
%             stimuli_size_local(l,:)  = stimuli_size_local(l,:)+( pi_dense_local(:,locations_trials(l,m)).*powers_trials(l))';
%         end
%        
%     end
% end
% evoked_cell = cell(n_trial,1);
% for i_trial = 1:n_trial
%     evoked_cell_index = 0; % 0: background evnets
%     for i_cell = 1:n_cell_temp
%         k = stimuli_size_local(i_trial, i_cell);
%         if k > k_minimum
%             evoked_cell_index = [evoked_cell_index i_cell];
%         end
%     end
%     evoked_cell{i_trial} = evoked_cell_index;
% end
% stimuli_size_temp = stimuli_size_local;
% %%
% cell_params.V_th = v_th_known*ones(n_cell_temp,1);
% cell_params.V_reset = v_reset_known*ones(n_cell_temp,1);
% 
% cell_params.gain = zeros(n_cell_temp,1);
% for i_cell = 1 : n_cell_temp
%     %cell_params.gain(i_cell) = l23_cells_for_sim(local_shape_gain(i_cell)).optical_gain;
%     cell_params.gain(i_cell) = mean_gain;
%     cell_params.gain_sd(i_cell)= std([l23_cells_for_sim.optical_gain]);
%      cell_params.g(i_cell) =  mean([l23_cells_for_sim.g]);
% end
% 
% 
% [event_rates,expected_all] =Intensity_v6(stimuli_size_temp,  n_stimuli_grid,n_grid_voltage,...
%     t_vect,t_factor,k_minimum,...
%     cell_params, funcs,...
%     I_stimuli,sd_range,delay_params,...
%     mpp);
%% Iterative updates
% gain_fits=zeros(n_cell_temp,maxit);
% gamma_fits=zeros(n_cell_temp,maxit);
% mu_fits=zeros(n_cell_temp,maxit);
% gain_old = 0.003*ones(n_cell_temp,1);
% gamma_old= 0.009*ones(n_cell_temp,1);
% mu_old = 2*ones(n_cell_temp,1);
% sigma_old = ones(n_cell_temp,1);
% 
% convergence_epsilon_outer=1e-3;
% normalized_change_outer=1;
% num_iter=1;
% soft_threshold =0.05;
% while (normalized_change_outer > convergence_epsilon_outer) & (num_iter < maxit)
%     num_iter = num_iter+1;
%     %----------------------------------------------%
%     % EM
%     
% [gamma_path mu_path sigma_path total_time soft_assignments bg_rate]= ...
%     EM_fullmodel_v3(mpp, ...
%     event_rates,...
%     evoked_cell,expected_all, ...
%     n_cell_temp, gamma_old, mu_old, sigma_old, ...
%     convergence_epsilon,f_background, mean_background, sigma_background, ...
%     sparsity, gamma_threshold,maxit,t_vect,use_size,background_update);
% 
% 
%     
%     % lif-glm updates:
%     % Use Monte-Carlo method to update lif-glm parameters based on soft
%     % assignments
%     % should be turned into a function
%     
%     % Reformat the soft assignments for each cell
%     
%     soft_assignments_by_cell = cell(n_cell_temp,1);
%     % record the trial, spike time, and soft assignments
%     for i_trial = 1:n_trial
%         n_event = length(mpp(i_trial).times);
%         cell_list = evoked_cell{i_trial};
%         if n_event >0
%             if length(cell_list)>1 % at least one cell beside the background
%                 for i_cell = 2:length(cell_list)
%                     if length(soft_assignments_by_cell{cell_list(i_cell)})==0
%                         %
%                         soft_assignments_by_cell{cell_list(i_cell)}=[i_trial*ones(1,n_event); mpp(i_trial).times;...
%                             soft_assignments{i_trial}(:,i_cell)'];
%                     else
%                         soft_assignments_by_cell{cell_list(i_cell)}=[soft_assignments_by_cell{cell_list(i_cell)}...
%                             [i_trial*ones(1,n_event); mpp(i_trial).times;soft_assignments{i_trial}(:,i_cell)']];
%                     end
%                 end
%                 
%             end
%         end
%     end
%     
%     
%     
%     num_MC_lifglm = 20;
%     
%     gains_sample = zeros(n_cell_temp,num_MC_lifglm);
%     gains_sd_sample= zeros(n_cell_temp,num_MC_lifglm);
%     for i_MC = 1:num_MC_lifglm
%         
%         % Consider an alternative way
%         % Choosing only the most reliable ones for each cell:
%         % Also, drawing only one event per cell
%         cell_data = cell(n_cell_temp,1);
%         for i_cell = 1:n_cell_temp
%             soft_temp =soft_assignments_by_cell{i_cell};
%             if length(soft_temp)<1
%                 cell_data{i_cell}=struct();
%                 cell_data{i_cell}.responses = [];
%                 cell_data{i_cell}.stims = [];
%                 
%             else
%                 trial_list=  unique(soft_temp(1,:));
%                 cell_data{i_cell}=struct();
%                 cell_data{i_cell}.responses = [];
%                 cell_data{i_cell}.stims = [];
%                 
%                 for i_trial = 1:length(trial_list)
%                     trial_idx= soft_temp(1,:)==trial_list(i_trial);
%                     soft_this_trial = soft_temp(:,trial_idx);
%                     
%                     prob_sum = sum(soft_this_trial(3,:));
%                     
%                     response_this_trial = zeros(1,length(I_e_vect));
%                     if prob_sum > soft_threshold
%                         if prob_sum > 1
%                             prob_sample = soft_this_trial(3,:)/prob_sum;
%                         else
%                             prob_sample=soft_this_trial(3,:);
%                         end
%                         
%                         prob_sample = [1- sum(prob_sample) prob_sample];
%                         r_temp = rand(1);
%                         i_event = min(find(r_temp<cumsum(prob_sample)));
%                         if i_event > 1
%                             response_this_trial(round(soft_this_trial(2,i_event-1)))=1;
%                         end
%                     end
%                     
%                     if sum(response_this_trial)>0
%                         cell_data{i_cell}.responses = [cell_data{i_cell}.responses; response_this_trial];
%                         cell_data{i_cell}.stims =  [cell_data{i_cell}.stims; ...
%                             I_e_vect*stimuli_size_temp(trial_list(i_trial), i_cell)];
%                     end
%                 end 
%             end
%         end
%         
%         lif_glm_gains= zeros(n_cell_temp,1);
%         for i_cell = 1:n_cell_temp
%             
%             N_cell = length(cell_data{i_cell}.responses);
%             in_params.g = cell_params.g(i_cell);
%             if (N_cell/length(I_e_vect))>0
%                 % LIF-GLM fits
%                 %-------------------------------------%
%                 responses=cell_data{i_cell}.responses;
%                 
%                 stims=cell_data{i_cell}.stims;
%                 
%                 [stats_conv] = fit_lifglm_v2(responses, ...
%                     stims,in_params,v_reset_known);
%                 %-------------------------------------%
%                 %lif_glm_gains(i_cell)=stats_conv.beta(2);
%                 lif_glm_gains(i_cell)=stats_conv.beta;
%                 
%             else
%                 lif_glm_gains(i_cell)=cell_params.gain(i_cell);
%                 gains_sample(i_cell,i_MC)=cell_params.gain_sd(i_cell);
%             end
%         end
%         
%         gains_sample(:,i_MC)=lif_glm_gains;
%     end
%     
%     % Evaluate the updates:
%     gamma_current= gamma_path(:,end);
%     sigma_current = sigma_path(:,end);
%     mu_current = mu_path(:,end);
%     gain_current = median(gains_sample,2);
%     for i_cell = 1:n_cell_temp
%         gain_sd_current(i_cell) = std(gains_sd_sample(i_cell,:));
%         if gain_sd_current(i_cell)==0
%             gain_sd_current(i_cell)=cell_params.gain_sd(i_cell);
%         end
%     end
%     
%     
%     
%     normalized_change_outer = norm(gamma_current - gamma_old)/(norm(gamma_old)+1) + norm(mu_current - mu_old)/(norm(mu_old)+1)+...
%         norm(sigma_current - sigma_old)/(norm(sigma_old)+1)+norm(gain_current-gain_old)/(norm(gain_old)+1);
%     
%     gamma_old =  gamma_current;
%     sigma_old = sigma_current;
%     mu_old = mu_current;
%     gain_old = gain_current;
%     
%     % Update the intensities
%     cell_params.V_th = v_th_known*ones(n_cell_temp,1);
%     cell_params.V_reset = v_reset_known*ones(n_cell_temp,1);
%     
%     % The local gains: update gains based on the lif-glm fits
%     cell_params.gain = gain_current;
%     cell_params.gain_sd = gain_sd_current;
%     
% 
% [event_rates,expected_all] =Intensity_v6(stimuli_size_temp,  n_stimuli_grid,n_grid_voltage,...
%     t_vect,t_factor,k_minimum,...
%     cell_params, funcs,...
%     I_stimuli,sd_range,delay_params,...
%     mpp);
%     
%    gain_fits(:,num_iter)=gain_current;
%    mu_fits(:,num_iter)=mu_current;
%    gamma_fits(:,num_iter)=gamma_current;
%    
%      fprintf('Changes %d\n', normalized_change_outer);
%   
%     
% end

% %%
% gamma_initial =zeros(n_cell_local,1);
% gamma_initial(stimulated_cells)=gamma_ini;
%  save('initial_fits_data_new.mat','gamma_initial');
% 
% %%
% gamma_all=zeros(n_cell_local,1);
% gamma_temp =zeros(sum(stimulated_cells),1);
% gamma_temp(connected_cells,1)=gamma_current;
% 
% gamma_all(stimulated_cells)=gamma_temp;
% %%
% 
% gain_all=zeros(n_cell_local,1);
% gain_temp =zeros(sum(stimulated_cells),1);
% gain_temp(connected_cells,1)=gain_current;
% 
% gain_all(stimulated_cells)=gain_temp;
%%
% save('final_fits_data_new.mat','gamma_all','gain_all');

%%
% %%
% figure(4)
% 
% temp2 = scatter(cell_locs(gamma_all>0.05,2),cell_locs(gamma_all>0.05,1),...
%     200*gamma_all(gamma_all>0.05));
% set(temp2,'MarkerFaceColor','r');
% alpha(temp2,0.5);
% 
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])
% 
% xlim([-150,150]);
% ylim([-150,150]);
% axis off;
% % outflnm =strcat('./');
% % saveas(4,strcat(outflnm,'Gamma_final_new','.jpg'));
% 
% 
% figure(5)
% temp2 = scatter(local_locations(gamma_all>0.05,2)+151,local_locations(gamma_all>0.05,1)+151,...
%     200*gamma_all(gamma_all>0.05));
% set(temp2,'MarkerFaceColor','r');
% alpha(temp2,0.5);
% hold on;
% 
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])
% 
% xlim([-13,311]);
% ylim([-13,311]);
% axis off;
% % Add legends to the plot 
% temp5 = scatter([0 0 0],[300 280 260],...
%     200*[1 0.5 0.2]);
% set(temp5,'MarkerFaceColor','r');
% alpha(temp5,0.5);
% 
% 
% txt4 = '\gamma = 1';
% txt5 = '\gamma = 0.5';
% txt6 = '\gamma = 0.2';
% 
% text(12,300,txt4)
% 
% text(12,280,txt5)
% 
% text(12,260,txt6)
% 
% rectangle('Position',[-13 250 72 60])
% 
% 
% outflnm =strcat('../Figures/Results/');
% saveas(5,strcat(outflnm,'Gamma_final_data','.jpg'));
