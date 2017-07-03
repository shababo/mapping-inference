%%
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Load real data
load('./Environments/6_3_s2c3_mpp_and_stim_data.mat')
Z_dense =target_locs;
powers_trials = stim_pow;
locations_trials=target_inds;
single_trial_limit=max(find(isnan(target_inds(:,2))));
locations_trials = target_inds(1:single_trial_limit,:);

powers_trials = stim_pow(1:single_trial_limit);
load('./Environments/6_3_cell_locs_reduced_and_fixed.mat')
local_locations = cell_locs;
n_trial = size(locations_trials,1);
n_cell_local=size(cell_locs,1);
%%
%------------------------------------------------------------------%
load('./Environments/chrome-template-3ms.mat');
downsamp=1;
max_time=300;
power_level = unique(powers_trials);
num_power_level=length(power_level);
I_stimuli=template(1:downsamp:max_time);
evoked_params.stim_start = 1;
evoked_params.stim_end = length(I_stimuli);
T=length(I_stimuli); % total time at 20k Hz
dt=1;
t_vect= dt:dt:T;
v_reset_known=-4e3;
min_time=80;
for i_trial = 1:n_trial
    if mpp(i_trial).num_events >0
        range_idx = mpp(i_trial).times<max_time & mpp(i_trial).times>min_time ;        
        mpp(i_trial).num_events = sum(range_idx);
        mpp(i_trial).times = mpp(i_trial).times(range_idx);
        mpp(i_trial).amp =mpp(i_trial).amp(range_idx);
    end
end

n_stimuli_grid=20;
n_grid_voltage=400;
t_factor=1;
gap_stimuli=0.5;
first_only=true;
stimulus_threshold=1e-3;
funcs.invlink = @invlink_test;%@(resp) log(1 + exp(resp));%@(x) exp(x);
n_grid_time = length(I_stimuli);

load('./Environments/l23_template_cell.mat');
%l23_average_shape
temp=l23_average_shape;
temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;
cell_params.locations =  local_locations;
cell_params.shape_gain = ones(n_cell_local,1);
shape_template = struct();
shape_template.shape= l23_average_shape;
[pi_dense_local, inner_normalized_products] = get_weights_v2(cell_params, shape_template,Z_dense);


stimuli_size_local=zeros(n_trial,n_cell_local);
for l = 1:n_trial
    for m = 1:size(locations_trials,2)
        if isnan(locations_trials(l,m))
        else
            stimuli_size_local(l,:) = stimuli_size_local(l,:)+(pi_dense_local(:,locations_trials(l,m)).*powers_trials(l))';
        end
    end
end
n_trial = size(stimuli_size_local,1);

%% Initialization
% Gamma:
% LIF-GLM parameters
%load('6_3_single_reduced.mat');
%%
%find(gamma_final_all>0.2)
max_radius = 10; % pick the maximum radius to consider 

%for this_cell = 1:n_cell_local % pick a cell
%
this_cell=71;
% Find stimuli locations that are close to this cell
% max_radius = 10; % pick the maximum radius to consider 
% events_prop = zeros(length(gamma_final_all),1);
% for this_cell = 1:length(gamma_final_all)
num_locs = size(target_locs,1);

diffs = target_locs - local_locations(this_cell,:);
relevant_locs = sqrt(sum(diffs.^2,2))<max_radius;
temp = 1:num_locs;
relevant_locs_index = temp(relevant_locs);
stim_counts = zeros(length(relevant_locs_index),1);
relevant_trials_index = [];
for i_trial = 1:n_trial
    if sum(relevant_locs_index == locations_trials(i_trial,1) )>0
        relevant_trials_index = [relevant_trials_index i_trial];
        stim_counts(find(relevant_locs_index == locations_trials(i_trial,1)))=...
            stim_counts(find(relevant_locs_index == locations_trials(i_trial,1)))+1;
    end
    if sum(relevant_locs_index == locations_trials(i_trial,2) )>0
        relevant_trials_index = [relevant_trials_index i_trial];
        stim_counts(find(relevant_locs_index == locations_trials(i_trial,1)))=...
            stim_counts(find(relevant_locs_index == locations_trials(i_trial,1)))+1;
    end
    
    if  sum(relevant_locs_index == locations_trials(i_trial,3) )>0
        relevant_trials_index = [relevant_trials_index i_trial];
        stim_counts(find(relevant_locs_index == locations_trials(i_trial,1)))=...
            stim_counts(find(relevant_locs_index == locations_trials(i_trial,1)))+1;
    end
    
end
if length(relevant_trials_index)>0
    
mpp_this_cell=mpp(relevant_trials_index);
locations_trials_this_cell = locations_trials(relevant_trials_index,:);
powers_trials_this_cell = powers_trials(relevant_trials_index);
stim_this_cell = stimuli_size_local(relevant_trials_index,this_cell);
% if size(mpp_this_cell,2) == 0
%     
% else
% events_prop(this_cell) = length([mpp_this_cell.times])/size(mpp_this_cell,2);
% end
% end
%stimuli_size=stimuli_size_local(relevant_trials_index,relevant_cells_index ) ;
%[~, this_cell]=max(events_prop);
% Draw the spike times 
figure(this_cell)
plot(sort([mpp_this_cell(powers_trials_this_cell==50).times]),1:length([mpp_this_cell(powers_trials_this_cell==50).times]),...
    '.','col','g','MarkerSize',20)
hold on;
plot(sort([mpp_this_cell(powers_trials_this_cell==75).times]),1:length([mpp_this_cell(powers_trials_this_cell==75).times]),...
    '.','col','b','MarkerSize',20)
hold on;
plot(sort([mpp_this_cell(powers_trials_this_cell==100).times]),1:length([mpp_this_cell(powers_trials_this_cell==100).times]),...
'.','col','r','MarkerSize',20)
hold off;
ylim([0,20]);
xlim([0,300]);
line([mean([mpp_this_cell(powers_trials_this_cell==50).times]) mean([mpp_this_cell(powers_trials_this_cell==50).times]) ], [0 20],'col','g')
line([mean([mpp_this_cell(powers_trials_this_cell==75).times]) mean([mpp_this_cell(powers_trials_this_cell==75).times]) ], [0 20],'col','b')
line([mean([mpp_this_cell(powers_trials_this_cell==100).times]) mean([mpp_this_cell(powers_trials_this_cell==100).times]) ], [0 20],'col','r')
xlabel('Time (1/20 ms)');
ylabel('Ranks (irrelevant)');

%saveas(this_cell,strcat('./sorted_time', num2str(this_cell),'.png'));
end
%end

%%
%%
%find(gamma_final_all>0.2)
max_radius = 10; % pick the maximum radius to consider 

%for this_cell = 1:n_cell_local % pick a cell
%
this_cell=71;
% Find stimuli locations that are close to this cell
% max_radius = 10; % pick the maximum radius to consider 
% events_prop = zeros(length(gamma_final_all),1);
% for this_cell = 1:length(gamma_final_all)
num_locs = size(target_locs,1);

diffs = target_locs - local_locations(this_cell,:);
relevant_locs = sqrt(sum(diffs.^2,2))<max_radius;
temp = 1:num_locs;
relevant_locs_index = temp(relevant_locs);
stim_counts = zeros(length(relevant_locs_index),1);
relevant_trials_index = [];
for i_trial = 1:n_trial
    if sum(relevant_locs_index == locations_trials(i_trial,1) )>0
        relevant_trials_index = [relevant_trials_index i_trial];
        stim_counts(find(relevant_locs_index == locations_trials(i_trial,1)))=...
            stim_counts(find(relevant_locs_index == locations_trials(i_trial,1)))+1;
    end
    if sum(relevant_locs_index == locations_trials(i_trial,2) )>0
        relevant_trials_index = [relevant_trials_index i_trial];
        stim_counts(find(relevant_locs_index == locations_trials(i_trial,1)))=...
            stim_counts(find(relevant_locs_index == locations_trials(i_trial,1)))+1;
    end
    
    if  sum(relevant_locs_index == locations_trials(i_trial,3) )>0
        relevant_trials_index = [relevant_trials_index i_trial];
        stim_counts(find(relevant_locs_index == locations_trials(i_trial,1)))=...
            stim_counts(find(relevant_locs_index == locations_trials(i_trial,1)))+1;
    end
    
end
    
mpp_this_cell=mpp(relevant_trials_index);
locations_trials_this_cell = locations_trials(relevant_trials_index,:);
powers_trials_this_cell = powers_trials(relevant_trials_index);
stim_this_cell = stimuli_size_local(relevant_trials_index,this_cell);
% if size(mpp_this_cell,2) == 0
%     
% else
% events_prop(this_cell) = length([mpp_this_cell.times])/size(mpp_this_cell,2);
% end
% end
%stimuli_size=stimuli_size_local(relevant_trials_index,relevant_cells_index ) ;
%[~, this_cell]=max(events_prop);
% Draw the spike times 
figure(this_cell)
plot(sort([mpp_this_cell(powers_trials_this_cell==50).times]),40+(1:length([mpp_this_cell(powers_trials_this_cell==50).times])),...
    '.','col','g','MarkerSize',20)

hold on;
plot(sort([mpp_this_cell(powers_trials_this_cell==75).times]),20+(1:length([mpp_this_cell(powers_trials_this_cell==75).times])),...
    '.','col','b','MarkerSize',20)

hold on;
plot(sort([mpp_this_cell(powers_trials_this_cell==100).times]),1:length([mpp_this_cell(powers_trials_this_cell==100).times]),...
'.','col','r','MarkerSize',20)
ylim([0,60]);
xlim([0,300]);

hold off;
xlabel('Time (1/20 ms)');
ylabel('Indices (irrelevant)');

saveas(this_cell,strcat('./sorted_time', num2str(this_cell),'.png'));

%%
% %%
% figure(1)
% % Locate  one isolated cell
% % this_cell = find(local_locations(:,2) < 50 & local_locations(:,2) > 30 & local_locations(:,1) < -50 & local_locations(:,2) > -80);
% % gamma_first_fit(this_cell)
% 
% temp1 = scatter(local_locations(:,2)+151,-local_locations(:,1)-151,...
%     40);
% set(temp1,'MarkerEdgeColor','g','MarkerFaceColor','g');
% alpha(temp1,1);
% hold on;
% 
% temp2 = scatter(local_locations(this_cell,2)+151,-local_locations(this_cell,1)-151,...
%     40,'MarkerEdgeColor','b','MarkerFaceColor','b',...
%     'MarkerFaceAlpha',1);
% hold on;
% % 
% % 
% % for i = 1:length(relevant_locs_index)
% %     if stim_counts(i)>0
% %         this_locs = relevant_locs_index(i);
% %         temp4 = scatter(target_locs(this_locs,2)+151,-target_locs(this_locs,1)-151,...
% %             stim_counts(i),'MarkerEdgeColor','r','MarkerFaceColor','r',...
% %             'MarkerFaceAlpha',0.5);
% %         hold on;
% %     end
% % end
% 
% %
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])
% 
% xlim([local_locations(this_cell,2)+151-20,local_locations(this_cell,2)+151+20]);
% ylim([-local_locations(this_cell,1)-151-20,-local_locations(this_cell,1)-151+20]);
% xlabel('x coordinate (micron)')
% ylabel('y coordinate (micron)')
% % saveas(1,'./Cell_71.png');
% 
% %% Use high power:
% power_idx =powers_trials_this_cell>50 ;% powers_trials_this_cell<100 ;
% 
% v_th_known=0;
% gain_initial =0;
% mean_delay = 35; %subtract the delays..
% in_params.g=mean([l23_cells_for_sim.g]);
% 
% responses = zeros(0,length(I_stimuli));
% responses_reg = zeros(0,length(I_stimuli));
% 
% stims = zeros(0,length(I_stimuli));
% counter = 0;
% for i = 1:length(mpp_this_cell)
%     if power_idx(i)
%         if length(mpp_this_cell(i).times)>0
%         counter = counter+1;
%         responses(counter,round(mpp_this_cell(i).times))=1;
%     times_adj=round(mpp_this_cell(i).times)-mean_delay;
%     responses_reg(counter,times_adj(times_adj>0))= 1;
%     stims(counter,:) = I_stimuli*stimuli_size_local(relevant_trials_index(i),this_cell);
%         end
%     end
% end
% 
% [stats_conv]=fit_lifglm_v4(responses_reg, stims,in_params,v_reset_known,v_th_known,first_only);
% if length(stats_conv.beta)>1
% gain_initial=stats_conv.beta(2)
% v_th_initial = -stats_conv.beta(1)
% else
% gain_initial=stats_conv.beta
%     v_th_initial = v_th_known
%     end
% %% Analysis:
% % 1, Initial fits
% %   - Estimate the firing rate given initial delay distribution and lif-glm
% %   parameters
% %   - Estimate the soft assignments and gammas given the fitted values
% 
% V_threshold = -50;
% cell_params.V_th = v_th_initial;
% cell_params.V_reset = v_reset_known;
% cell_params.gain = gain_initial;
% cell_params.gain_sd= std([l23_cells_for_sim.optical_gain]);
% cell_params.g =  in_params.g;
% 
% n_delay_grid = 200;
% outputM=false;
% delay_params_est.type=1;
% delay_params_est.mean=mean_delay;
% delay_params_est.std=20;
% 
% stimulus_threshold=1e-3;
% [Stimuli_grid, Intensity_grid]=Intensity_v8(...
%     stimuli_size_local(relevant_trials_index,this_cell), mpp_this_cell,I_stimuli,... % data from exp
%     cell_params,... % estimated parameters
%     funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
%     n_stimuli_grid,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
%     V_threshold,stimulus_threshold,first_only);
% %%
% for i=1:3
% figure(i)
% plot(Intensity_grid{i})
% end
% 
% %%
% sd_range=0.05;
% n_trial = length(mpp_this_cell);
% est_intensity = cell(length(relevant_trials_index),1);
% %-----------------------------------------%
% % Convolute the intensity with delay distribution
% delay_prob = zeros(2*n_delay_grid+1,1);
% if delay_params_est.std == 0
%     delay_prob(n_delay_grid+1)=1;
%     min_delay=0;
%     max_delay=0;
% else
%     delay_prob = normpdf( -(-n_delay_grid:n_delay_grid),delay_params_est.mean,...
%         delay_params_est.std);
%     % we approximate the probability with densities
%     delay_prob = delay_prob/sum(delay_prob);
%     min_delay = 1-1-n_delay_grid;
%     max_delay = length(delay_prob)-1-n_delay_grid;
% end
% 
% for i_stimuli = 1:length(Stimuli_grid)
%     M_grid_intensity{i_stimuli}=zeros(length(Intensity_grid{1}),1);
%     for i_t = 1:length(Intensity_grid{1})
%         idx_time = max(i_t+min_delay,1): min(i_t+max_delay,n_grid_time);
%         idx_delay = (min(idx_time)-i_t+n_delay_grid+1) : (max(idx_time)-i_t+n_delay_grid+1);
%         M_grid_intensity{i_stimuli}(i_t)=delay_prob(idx_delay)*Intensity_grid{i_stimuli}(idx_time);
%     end
% end
% 
% %%
% power_levels =[50 75 100];
% for i = 1:3
%     i_power = power_levels(i);
% idx=powers_trials_this_cell==i_power;   
% temp =  M_grid_intensity{i};
% 
% figure(i)
% 
%     plot( (temp+background_rate) *10*sum(idx),'r','LineWidth',4);
%     hold on;
% counts =[mpp_this_cell(idx).times];    
% %freq = counts/300/sum(powers_trials_this_cell==i_power);
% histogram(counts,'BinWidth',10)
% % plot(freq,'col','b','LineWidth',4);
% length(counts)
% sum(temp+background_rate)*sum(idx)
% ylim([0,6]);
% hold off
% end
% 
% 
% %% Iteratively update the delay variable and delay distribution:
% 
% V_threshold = -50;
% cell_params.V_th = v_th_initial;
% cell_params.V_reset = v_reset_known;
% cell_params.gain = gain_initial;
% cell_params.gain_sd= std([l23_cells_for_sim.optical_gain]);
% cell_params.g =  in_params.g;
% 
% n_delay_grid = 200;
% outputM=false;
% gain_fit = gain_initial;
% v_th_fit= v_th_initial;
% %%
% 
% delay_params_est.type=1;
% delay_params_est.mean=mean_delay;
% delay_params_est.std=10;
% 
% epsilon_threshold = 1e-2;
% normalized_change = 1;
% while (normalized_change > epsilon_threshold) 
% cell_params.gain = gain_fit;
% cell_params.V_th = v_th_fit;
% [Stimuli_grid, Intensity_grid]=Intensity_v8(...
%     stimuli_size_local(relevant_trials_index,this_cell), mpp_this_cell,I_stimuli,... % data from exp
%     cell_params,... % estimated parameters
%     funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
%     n_stimuli_grid,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
%     V_threshold,stimulus_threshold,first_only);
% 
% 
% 
% % responses and stims have been defined
% responses_reg=responses;responses_reg(:,:)=0;
% delays=[];
% for i_trial = 1:size(stims,1)
%     k_temp = max(stims(i_trial,:))/max(I_stimuli)*gain_fit;
%     spikes=find(responses(i_trial,:));
%     
%     [~, idx_min]=min( abs(k_temp-Stimuli_grid));
%     intensity_temp=Intensity_grid{idx_min};
%     if length(spikes)>0
%         spike_first = spikes(1);
%         delay_seq =spike_first-(1:length(I_stimuli));
%         if delay_params_est.std > 0.5
%         else 
%             delay_params_est.std = 0.5;
%         end  
%         delay_prob = normpdf(delay_seq,...
%             delay_params_est.mean,delay_params_est.std);
%             
%         prod_prob = delay_prob.*intensity_temp';
%         [~, spike_max]= max(prod_prob);
%         responses_reg(i_trial,spike_max)=1;
%         delays =[delays spike_first-spike_max];
%     end
%     
% end
% 
% [stats_conv]=fit_lifglm_v4(responses_reg, stims,in_params,v_reset_known,v_th_known,first_only);
% 
% if v_th_known == 0
% normalized_change = abs(stats_conv.beta(2)-gain_fit) /gain_fit + ...
%     abs(-stats_conv.beta(1)-v_th_fit) /v_th_fit;%+...
%    % abs(delay_params_est.mean-mean(delays))/mean(delays);
% fprintf('Changes %d\n', normalized_change);
%     
% gain_fit=stats_conv.beta(2);
% v_th_fit= - stats_conv.beta(1);
% else 
% normalized_change =    abs(stats_conv.beta-gain_fit) /gain_fit+...
%    abs(delay_params_est.mean-mean(delays))/mean(delays);
% fprintf('Changes %d\n', normalized_change);
%     
% gain_fit=stats_conv.beta;
% end    
% % delay_params_est.mean=mean(delays);
% % delay_params_est.std=std(delays);
% end
% 
% %% Summary
% 
% sd_range=0.05;
% n_trial = length(mpp_this_cell);
% est_intensity = cell(length(relevant_trials_index),1);
% %-----------------------------------------%
% % Convolute the intensity with delay distribution
% delay_prob = zeros(2*n_delay_grid+1,1);
% if delay_params_est.std == 0
%     delay_prob(n_delay_grid+1)=1;
%     min_delay=0;
%     max_delay=0;
% else
%     delay_prob = normpdf( -(-n_delay_grid:n_delay_grid),delay_params_est.mean,...
%         delay_params_est.std);
%     % we approximate the probability with densities
%     delay_prob = delay_prob/sum(delay_prob);
%     min_delay = 1-1-n_delay_grid;
%     max_delay = length(delay_prob)-1-n_delay_grid;
% end
% 
% for i_stimuli = 1:length(Stimuli_grid)
%     M_grid_intensity{i_stimuli}=zeros(length(Intensity_grid{1}),1);
%     for i_t = 1:length(Intensity_grid{1})
%         idx_time = max(i_t+min_delay,1): min(i_t+max_delay,n_grid_time);
%         idx_delay = (min(idx_time)-i_t+n_delay_grid+1) : (max(idx_time)-i_t+n_delay_grid+1);
%         M_grid_intensity{i_stimuli}(i_t)=delay_prob(idx_delay)*Intensity_grid{i_stimuli}(idx_time);
%     end
% end
% 
% 
% %%
% figure(4)
% histogram(delays)
% 
% 
% power_levels =[50 75 100];
% for i = 1:3
% i_power = power_levels(i);
% idx=powers_trials_this_cell==i_power;   
% temp =  M_grid_intensity{i};
% 
% figure(i)
% 
%     plot( (temp+background_rate) *10*sum(idx),'r','LineWidth',4);
%     hold on;
% counts =[mpp_this_cell(idx).times];    
% %freq = counts/300/sum(powers_trials_this_cell==i_power);
% histogram(counts,'BinWidth',10)
% % plot(freq,'col','b','LineWidth',4);
% length(counts)
% sum(temp+background_rate)*sum(idx)
% ylim([0,6]);
% hold off
% end
