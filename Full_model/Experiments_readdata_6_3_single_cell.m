%%
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% load real data
load('./Environments/6_3_s2c3_mpp_and_stim_data.mat')
%target_locs,
%% Setting experiment details based on the data
Z_dense =target_locs;
powers_trials = stim_pow;
locations_trials=target_inds;
single_trial_limit=max(find(isnan(target_inds(:,2))));
locations_trials = target_inds(1:single_trial_limit,:);

powers_trials = stim_pow(1:single_trial_limit);
load('./Environments/6_3_cell_locs_reduced_and_fixed.mat')
local_locations = cell_locs;
n_trial = size(locations_trials,1);
%%
%------------------------------------------------------------------%
trials_specific_variance= 0;
load('./Environments/l23_template_cell.mat');
temp=l23_average_shape;
temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;
load('./Environments/l23_cells_for_sim.mat');
mean_gain = mean([l23_cells_for_sim.optical_gain])/2;
n_trial = size(locations_trials,1);
n_cell_local = size(local_locations,1);
cell_params.locations =  local_locations;
cell_params.shape_gain = ones(n_cell_local,1);
shape_template = struct();
shape_template.shape= l23_average_shape;
[pi_dense_local, inner_normalized_products] = get_weights_v2(cell_params, shape_template,Z_dense);
load('./Environments/chrome-template-3ms.mat');
downsamp=1;
max_time=300;
power_level = unique(powers_trials);
num_power_level=length(power_level);
current_template=template(1:downsamp:max_time);
I_stimuli =current_template;
evoked_params.stim_start = 1;
evoked_params.stim_end = length(I_stimuli);
data_params.T = length(I_stimuli); % total time at 20k Hz
data_params.dt = 1; % 1/20 ms
k_minimum = 0.1; % minimum stimulus intensity to consider
v_reset_known=-4e3;
cell_params.V_th = 15*ones([n_cell_local,1]);
cell_params.V_reset = v_reset_known*ones([n_cell_local,1]);
max_time = 300;
min_time=80;
for i_trial = 1:n_trial
    if mpp(i_trial).num_events >0
        range_idx = mpp(i_trial).times<max_time & mpp(i_trial).times>min_time ;
        
        mpp(i_trial).num_events = sum(range_idx);
        mpp(i_trial).times = mpp(i_trial).times(range_idx);
        mpp(i_trial).amp =mpp(i_trial).amp(range_idx);
    end
    
end
sigma_unknown=1;
T=length(I_stimuli); % total time at 20k Hz
dt=1;
t_vect= dt:dt:T;

n_stimuli_grid=20;
n_grid_voltage=400;
t_factor=1;
gap_stimuli=0.5;
first_only=true;
stimulus_threshold=0.1;
funcs.invlink = @invlink_test;%@(resp) log(1 + exp(resp));%@(x) exp(x);
n_grid_time = length(I_stimuli);
%
stimuli_size_local=zeros(n_trial,n_cell_local);
for l = 1:n_trial
    for m = 1:size(locations_trials,2)
        if isnan(locations_trials(l,m))
        else
            stimuli_size_local(l,:) = stimuli_size_local(l,:)+( pi_dense_local(:,locations_trials(l,m)).*powers_trials(l))';
        end
    end
end
n_trial = size(stimuli_size_local,1);
%----------------------------------------------------------%
%%
load('6_3_single_reduced.mat');
%% Pick the cell and relevant trials:
%---------------------------%
find(gamma_final_all>0.2)
this_cell=30; % pick a cell
max_radius = 40; % pick the maximum radius to consider 
%%
% Find stimuli locations that are close to this cell

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

%%
figure(1)
% Locate  one isolated cell
% this_cell = find(local_locations(:,2) < 50 & local_locations(:,2) > 30 & local_locations(:,1) < -50 & local_locations(:,2) > -80);
% gamma_first_fit(this_cell)
temp2 = scatter(local_locations(this_cell,2)+151,-local_locations(this_cell,1)-151,...
    100,'MarkerEdgeColor','g','MarkerFaceColor','b',...
    'MarkerFaceAlpha',0.5);
hold on;

temp1 = scatter(local_locations(:,2)+151,-local_locations(:,1)-151,...
    2);
set(temp1,'MarkerEdgeColor','g','MarkerFaceColor','g');
alpha(temp1,1);
hold on;


for i = 1:length(relevant_locs_index)
    if stim_counts(i)>0
        this_locs = relevant_locs_index(i);
        temp4 = scatter(target_locs(this_locs,2)+151,-target_locs(this_locs,1)-151,...
            stim_counts(i),'MarkerEdgeColor','r','MarkerFaceColor','r',...
            'MarkerFaceAlpha',0.5);
        hold on;
    end
end

%
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])

xlim([-20,313]);
ylim([-313,20]);


%% Pick related cells 
relevant_cells_index = [];
stim_threshold = 10;
for i_cell = 1:n_cell_local
    
 if sum(stimuli_size_local(relevant_trials_index,i_cell)>stim_threshold)>0
    relevant_cells_index=[ relevant_cells_index i_cell]; 
    
 end
end

%% Only consider the cell of interest
mpp_this_cell=mpp(relevant_trials_index);
locations_trials_this_cell = locations_trials(relevant_trials_index,:);
powers_trials_this_cell = powers_trials(relevant_trials_index);
stimuli_size=stimuli_size_local(relevant_trials_index,relevant_cells_index ) ;
%%---------------------------------------------%%

%%

V_threshold = -50;
cell_params.V_th = 15*ones(length(relevant_cells_index),1);
cell_params.V_reset = v_reset_known*ones(length(relevant_cells_index),1);
cell_params.gain = gain_final_all(relevant_cells_index);
cell_params.gain_sd= std([l23_cells_for_sim.optical_gain])*ones(length(relevant_cells_index),1);
cell_params.g =  mean([l23_cells_for_sim.g])*ones(length(relevant_cells_index),1);

n_delay_grid = 200;
outputM=false;
delay_params_est.type=1;
delay_params_est.mean=35;
delay_params_est.std=10;
sd_range=0.5;

[Stimuli_grid, Intensity_grid]=Intensity_v8(...
    stimuli_size, mpp_this_cell,I_stimuli,... % data from exp
    cell_params,... % estimated parameters
    funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
    n_stimuli_grid,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
    V_threshold,stimulus_threshold,first_only);


%%
n_trial = length(mpp_this_cell);
intensities_relevant_cells = cell(n_trial,length(relevant_cells_index));

delay_prob = zeros(2*n_delay_grid+1,1);
if delay_params_est.std == 0
    delay_prob(n_delay_grid+1)=1;
    min_delay=0;
    max_delay=0;
else
    delay_prob = normpdf( -(-n_delay_grid:n_delay_grid),delay_params_est.mean,...
        delay_params_est.std);
    % we approximate the probability with densities
    delay_prob = delay_prob/sum(delay_prob);
    min_delay = 1-1-n_delay_grid;
    max_delay = length(delay_prob)-1-n_delay_grid;
end

for i = 1:length(relevant_cells_index)
    i_cell = relevant_cells_index(i);
for i_stimuli = 1:length(Stimuli_grid)
    M_grid_intensity{i_stimuli}=zeros(length(Intensity_grid{1}),1);
    for i_t = 1:length(Intensity_grid{1})
        idx_time = max(i_t+min_delay,1): min(i_t+max_delay,n_grid_time);
        idx_delay = (min(idx_time)-i_t+n_delay_grid+1) : (max(idx_time)-i_t+n_delay_grid+1);
        M_grid_intensity{i_stimuli}(i_t)=delay_prob(idx_delay)*Intensity_grid{i_stimuli}(idx_time);
    end
end
%------------------------------------------------%
for i_trial = 1:n_trial
    %------------------------------------------------%
    % Allowing for noise in the gain estimates
    k_up = stimuli_size(i_trial,i)*(cell_params.gain(i)+sd_range*cell_params.gain_sd(i));
    k_low =stimuli_size(i_trial,i)*(cell_params.gain(i)-sd_range*cell_params.gain_sd(i));
    intensity_temp = zeros([length(t_vect) 1]);
    index_seq = (k_low<Stimuli_grid &  k_up>Stimuli_grid);
    if sum(index_seq)>0
        for i_grid = 1:length(Stimuli_grid)
            if index_seq(i_grid)>0
                intensity_temp= intensity_temp+M_grid_intensity{i_grid};
            end
        end
        intensity_temp=intensity_temp/sum(index_seq);
        intensities_relevant_cells{i_trial,i}=intensity_temp;
    end
    %------------------------------------------------%
end
end

%% Plot the intensities against the observed values 
i_trial =200;
if length(mpp_this_cell(i_trial).times)>0

plot(mpp_this_cell(i_trial).times,zeros(1,length(mpp_this_cell(i_trial).times)),'.','col','b',...
    'MarkerSize',20);
hold on;
for i=1:length(relevant_cells_index)
    if relevant_cells_index(i)==this_cell
%         plot(intensities_relevant_cells{i_trial,i},'-','col','r'); %
%         hold on; %         intensity only

        plot(intensities_relevant_cells{i_trial,i}*gamma_final_all(relevant_cells_index(i)),'-','col','r');
        hold on;% with gammas
    else
%         plot(intensities_relevant_cells{i_trial,i},'--','col','k');
%         hold on;
        plot(intensities_relevant_cells{i_trial,i}*gamma_final_all(relevant_cells_index(i)),'--','col','k');
        hold on;
    end
end
hold off;
ylim([0 max(max([intensities_relevant_cells{i_trial,:}]))] );

% The soft assignments for the events:
soft_assignments_labeled{relevant_trials_index(i_trial)}(:,this_cell+1)
else
    fprintf('No events in this trial.\n');
end
%%
% 
%% Initialization
% Gamma:
% LIF-GLM parameters
find(gamma_final_all>0.2)
this_cell=94; % pick a cell
max_radius = 20; % pick the maximum radius to consider 
%%
% Find stimuli locations that are close to this cell

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
%%
relevant_trials_index
mpp_this_cell=mpp(relevant_trials_index);
locations_trials_this_cell = locations_trials(relevant_trials_index,:);
powers_trials_this_cell = powers_trials(relevant_trials_index);
%stimuli_size=stimuli_size_local(relevant_trials_index,relevant_cells_index ) ;

%%
gain_initial =(0);
mean_delay = 35; %subtract the delays..
in_params.g=mean([l23_cells_for_sim.g]);

responses = zeros(length(mpp_this_cell),length(I_stimuli));
stims = zeros(length(mpp_this_cell),length(I_stimuli));
for i = 1:length(mpp_this_cell)
    times_adj=round(mpp_this_cell(i).times)-mean_delay;
    responses(i,times_adj(times_adj>0))= 1;
    stims(i,:) = I_stimuli*stimuli_size_local(relevant_trials_index(i),this_cell);
end

[stats_conv]=fit_lifglm_v3(responses, stims,in_params,v_reset_known,first_only);
gain_initial=stats_conv.beta;

%% Analysis:
% 1, Initial fits
%   - Estimate the firing rate given initial delay distribution and lif-glm
%   parameters
%   - Estimate the soft assignments and gammas given the fitted values

V_threshold = -50;
cell_params.V_th = 15;
cell_params.V_reset = v_reset_known;
cell_params.gain = gain_initial;
cell_params.gain_sd= std([l23_cells_for_sim.optical_gain]);
cell_params.g =  in_params.g;

n_delay_grid = 200;
outputM=false;
delay_params_est.type=1;
delay_params_est.mean=35;
delay_params_est.std=10;


[Stimuli_grid, Intensity_grid]=Intensity_v8(...
    stimuli_size_local(relevant_trials_index,this_cell), mpp_this_cell,I_stimuli,... % data from exp
    cell_params,... % estimated parameters
    funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
    n_stimuli_grid,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
    V_threshold,stimulus_threshold,first_only);
%%
for i=1:3
figure(i)
plot(Intensity_grid{i})
end

sd_range=0.2;
n_trial = length(mpp_this_cell);
est_intensity = cell(length(relevant_trials_index),1);
%-----------------------------------------%
% Convolute the intensity with delay distribution
delay_prob = zeros(2*n_delay_grid+1,1);
if delay_params_est.std == 0
    delay_prob(n_delay_grid+1)=1;
    min_delay=0;
    max_delay=0;
else
    delay_prob = normpdf( -(-n_delay_grid:n_delay_grid),delay_params_est.mean,...
        delay_params_est.std);
    % we approximate the probability with densities
    delay_prob = delay_prob/sum(delay_prob);
    min_delay = 1-1-n_delay_grid;
    max_delay = length(delay_prob)-1-n_delay_grid;
end

for i_stimuli = 1:length(Stimuli_grid)
    M_grid_intensity{i_stimuli}=zeros(length(Intensity_grid{1}),1);
    for i_t = 1:length(Intensity_grid{1})
        idx_time = max(i_t+min_delay,1): min(i_t+max_delay,n_grid_time);
        idx_delay = (min(idx_time)-i_t+n_delay_grid+1) : (max(idx_time)-i_t+n_delay_grid+1);
        M_grid_intensity{i_stimuli}(i_t)=delay_prob(idx_delay)*Intensity_grid{i_stimuli}(idx_time);
    end
end

%------------------------------------------------%
for i_trial = 1:n_trial
    %------------------------------------------------%
    % Allowing for noise in the gain estimates
    k_up = stimuli_size_local(relevant_trials_index(i_trial),this_cell)*(cell_params.gain+sd_range*cell_params.gain_sd);
    k_low =stimuli_size_local(relevant_trials_index(i_trial),this_cell)*(cell_params.gain-sd_range*cell_params.gain_sd);
    intensity_temp = zeros([length(t_vect) 1]);
    index_seq = (k_low<Stimuli_grid &  k_up>Stimuli_grid);
    if sum(index_seq)>0
        for i_grid = 1:length(Stimuli_grid)
            if index_seq(i_grid)>0
                intensity_temp= intensity_temp+M_grid_intensity{i_grid};
            end
        end
        intensity_temp=intensity_temp/sum(index_seq);
        expected_all(i_trial)=sum(intensity_temp(min_time:end));
        if length(mpp(i_trial).times)>0
            event_rates{i_trial}=intensity_temp( max(1,round(mpp_this_cell(i_trial).times)) );
        end
    end
    est_intensity{i_trial} =intensity_temp;
    %------------------------------------------------%
end

%%
power_levels =[50 75 100];
for i = 1:3
    i_power = power_levels(i);
idx=powers_trials_this_cell==i_power;   
temp = est_intensity{idx};

figure(i)

    plot( (temp{1}+background_rate) *10*sum(idx),'r','LineWidth',4);
    hold on;
counts =[mpp_this_cell(idx).times];    
%freq = counts/300/sum(powers_trials_this_cell==i_power);
histogram(counts,'BinWidth',10)
% plot(freq,'col','b','LineWidth',4);
hold off
end