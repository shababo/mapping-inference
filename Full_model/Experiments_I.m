% Loading functions and Data generation
% cd('C:/Users/Shizhe/Documents/GitHub/mapping-inference/Full_model/')
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Paramters in the simulations
trials_specific_variance= 0;
%% Loading templates from real data
load('./Environments/l23_cells_for_sim.mat');
l23_cells_for_sim_origin=l23_cells_for_sim;

% normalized the cell shapes
counter=0;

for i = 1:length(l23_cells_for_sim_origin)
   % if l23_cells_for_sim_origin(i).optical_gain>0.015 & l23_cells_for_sim_origin(i).optical_gain<0.02
        counter=counter+1;
        if counter==1
        templates=l23_cells_for_sim_origin(i);
        else
        templates(counter)=l23_cells_for_sim_origin(i);
        end
        temp=l23_cells_for_sim(i).shape;
        temp_max = max(max(max(temp)));
        templates(counter).shape = temp/temp_max;
       % templates(counter).optical_gain=0.02;%templates(counter).optical_gain;
    %end
end
l23_cells_for_sim=templates;
num_types_cell = length(l23_cells_for_sim);

load('./Environments/l23_template_cell.mat');
%l23_average_shape
temp=l23_average_shape;
temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;
%% load real data
load('./Environments/6_3_s2c3_mpp_and_stim_data.mat')
%% Setting experiment details based on the data
local_locations= cell_locs;
Z_dense =target_locs;
powers_trials = stim_pow;
locations_trials=target_inds;
%% Using only the single spot data 
single_trial_limit=max(find(isnan(target_inds(:,2))));
locations_trials = target_inds(1:single_trial_limit,:);
powers_trials = stim_pow(1:single_trial_limit);

%% Set seed for reproducibility
rng(12242,'twister');
n_trial = size(locations_trials,1);
% n_trial = 8000;
n_cell_local = size(local_locations,1);

local_shape_gain=randsample(1:num_types_cell,n_cell_local,true);
%% Use the estimates from the real data:
load('final_fits_6_3_single.mat');
local_gamma = gamma_all;
local_gamma(gamma_all<0.05)=0;
local_gamma(gamma_all>0.05 &gamma_all<0.2)=0.2;
local_connected = local_gamma>0;
local_gain=gain_all;
local_gain(208)=0.02;
local_gain(211)=0.015;
local_gain(325)=0.02;

%gain_all
%% Calculate the background rate using the low power stimulations
mpp_copy=mpp;
histogram([mpp_copy.times]);
% Why are there a lot of events around 500?
trial_counts = 0;
event_counts = 0;
for i_trial = 1:n_trial
    %if stim_pow(i_trial)==50
    event_counts = event_counts +sum(mpp_copy(i_trial).times>20 & mpp_copy(i_trial).times<80);
    %end
end
background_rate = event_counts/n_trial/60;
%%
cell_params.locations = local_locations;
cell_params.shape_gain = local_shape_gain';
cell_params.gains=local_gain;
shape_template= l23_cells_for_sim;
[pi_dense_all,~] = get_weights_v2(cell_params, shape_template,Z_dense);
%%
cell_params.locations =  local_locations;
cell_params.shape_gain = ones(n_cell_local,1);
shape_template = struct();
shape_template.shape= l23_average_shape;
[pi_dense_local, inner_normalized_products] = get_weights_v2(cell_params, shape_template,Z_dense);
%% Loading the current template using new template
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

maxit=100;
convergence_epsilon=1e-3;
mean_background = 2;
sigma_background  = 1;
%gamma_threshold=0.1;
k_minimum = 0.1; % minimum stimulus intensity to consider
%%
v_reset_known=-4e3;

cell_params.V_th = 15*ones([n_cell_local,1]);
cell_params.V_reset = v_reset_known*ones([n_cell_local,1]);
cell_params.gamma = local_gamma;
cell_params.amplitudes = abs(normrnd(3,1,[n_cell_local, 1]));
cell_params.sigma_across =  abs(normrnd(0.5,0.1,[n_cell_local, 1]));
cell_params.sigma_within = abs(normrnd(0.5,0.1,[n_cell_local, 1]));
cell_params.locations = local_locations;
cell_params.shape_gain = local_shape_gain;
cell_params.gains=local_gain;

bg_params.mean = 2;
bg_params.sigma = 1.5;
bg_params.firing_rate = background_rate; %spike/sec
shape_template= l23_cells_for_sim;

delay_params.type=1;
delay_params.mean=35;
delay_params.std=15;
%%
rng(12242,'twister');
[mpp_temp, ~, ~] = generate_data_v3(...
    locations_trials(1:n_trial,:),powers_trials(1:n_trial,:),...
    pi_dense_all,k_minimum,cell_params, shape_template, ...
    I_stimuli, data_params,bg_params,trials_specific_variance,...
    delay_params);
mpp=mpp_temp;
%%
%load('mpp_sim_6_3_single.mat','mpp');
    max_time = 300;
min_time=60;
for i_trial = 1:n_trial
    if mpp_copy(i_trial).num_events >0
        range_idx = mpp_copy(i_trial).times<max_time & mpp_copy(i_trial).times>min_time ;

        mpp_copy(i_trial).num_events = sum(range_idx);
        mpp_copy(i_trial).times = mpp_copy(i_trial).times(range_idx);
        mpp_copy(i_trial).amp =mpp_copy(i_trial).amp(range_idx);
    end

end
%
%     %% count the number of events:
   sum([mpp.assignments]==0)
   300*n_trial*background_rate
length([mpp.assignments])
length([mpp_copy(1:n_trial).times])
save('mpp_sim_6_3_single.mat','mpp');

%%
stimuli_size_local=zeros(n_trial,n_cell_local);
for l = 1:n_trial
    for m = 1:size(locations_trials,2)
        if isnan(locations_trials(l,m))
        else
            stimuli_size_local(l,:) = stimuli_size_local(l,:)+( pi_dense_local(:,locations_trials(l,m)).*powers_trials(l))';
        end
    end
end
%% Check the qualities of the simulations:
connected_list = find(local_gamma>0);
all_assignments =[mpp.assignments];
for i_cell_idx = 1:length(connected_list)
    i_cell=connected_list(i_cell_idx);
    [i_cell sum(all_assignments==i_cell) sum(stimuli_size_local(:,i_cell)>20)]
end
%%
%------------------------------------%
% Estimating the marginal firing rate
sigma_unknown=1;
dt=1;
t_vect= dt:dt:max_time;
sd_range=1.5;
n_stimuli_grid=20;
n_grid_voltage=400;
t_factor=1;
gap_stimuli=0.5;
first_only=true;
stimulus_threshold=0.1;
funcs.invlink = @invlink_test;%@(resp) log(1 + exp(resp));%@(x) exp(x);
n_grid_time = length(I_stimuli);

%% Throw away cells that are not stimulated enough:
stim_threshold = 30;
stimulated_cells = sum(stimuli_size_local>stim_threshold )>5;
%%
stimuli_size_stimulated = stimuli_size_local(:,stimulated_cells);
n_cell_stimulated = sum(stimulated_cells);
evoked_cell_stimulated = cell(n_trial,1);
for i_trial = 1:n_trial
    evoked_cell_index = 0; % 0: background evnets
    for i_cell = 1:n_cell_stimulated
        k = stimuli_size_stimulated(i_trial, i_cell);
        if k > k_minimum
            evoked_cell_index = [evoked_cell_index i_cell];
        end
    end
    evoked_cell_stimulated{i_trial} = evoked_cell_index;
end
%% Reformat the mpp:
for i = 1:n_trial
    mpp(i).amp=mpp(i).amplitudes;
    mpp(i).times=mpp(i).event_times;
end
%% Analysis:
% 1, Initial fits
%   - Estimate the firing rate given initial delay distribution and lif-glm
%   parameters
%   - Estimate the soft assignments and gammas given the fitted values
V_threshold = -50;
cell_params.V_th = 15*ones(sum(stimulated_cells),1);
cell_params.V_reset = v_reset_known*ones(sum(stimulated_cells),1);
cell_params.gain = ones(n_cell_stimulated,1)*mean([l23_cells_for_sim.optical_gain]);
cell_params.gain_sd= ones(n_cell_stimulated,1)*std([l23_cells_for_sim.optical_gain]);
cell_params.g =  ones(n_cell_stimulated,1)*mean([l23_cells_for_sim.g]);

n_delay_grid = 200;
outputM=false;
delay_params_est.type=1;
delay_params_est.mean=35*ones(n_cell_stimulated,1);
delay_params_est.std=15*ones(n_cell_stimulated,1);


[Stimuli_grid, Intensity_grid]=Intensity_v8(...
    stimuli_size_stimulated, mpp,I_stimuli,... % data from exp
    cell_params,... % estimated parameters
    funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
    n_stimuli_grid,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
    V_threshold,stimulus_threshold,first_only);

expected_all = zeros(n_trial,n_cell_stimulated);
event_rates = cell(n_trial,1);
for i_trial = 1:n_trial
    if length(mpp(i_trial).times)>0
        event_rates{i_trial}=zeros(length(mpp(i_trial).times),n_cell_stimulated);
    end
end
for i_cell = 1:n_cell_stimulated
    
    %-----------------------------------------%
    % Convolute the intensity with delay distribution
    delay_prob = zeros(2*n_delay_grid+1,1);
    if delay_params_est.std(i_cell) == 0
        delay_prob(n_delay_grid+1)=1;
        min_delay=0;
        max_delay=0;
    else
        delay_prob = normpdf( -(-n_delay_grid:n_delay_grid),delay_params_est.mean(i_cell),delay_params_est.std(i_cell));
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
        k_up = stimuli_size_stimulated(i_trial, i_cell)*(cell_params.gain(i_cell)+sd_range*cell_params.gain_sd(i_cell));
        k_low = stimuli_size_stimulated(i_trial, i_cell)*(cell_params.gain(i_cell)-sd_range*cell_params.gain_sd(i_cell));
        intensity_temp = zeros([length(t_vect) 1]);
        index_seq = (k_low<Stimuli_grid &  k_up>Stimuli_grid);
        if sum(index_seq)>0
            for i_grid = 1:length(Stimuli_grid)
                if index_seq(i_grid)>0
                    intensity_temp= intensity_temp+M_grid_intensity{i_grid};
                end
            end
            intensity_temp=intensity_temp/sum(index_seq);
            expected_all(i_trial,i_cell)=sum(intensity_temp);
            if length(mpp(i_trial).times)>0
                event_rates{i_trial}(:,i_cell)=intensity_temp( max(1,round(mpp(i_trial).times)) );
            end
        end
        %------------------------------------------------%
    end
    fprintf('%d\n',i_cell);
end
%%
background_update=0;
f_background = background_rate;

mean_background = 0; %not actually used
sigma_background  =1;%not actually used
gamma_old= 0.009*ones(n_cell_stimulated,1);
mu_old = 2*ones(n_cell_stimulated,1);
sigma_old = ones(n_cell_stimulated,1);

use_size =0;
convergence_epsilon = 0.01;
maxit = 100;

sparsity_params.threshold=4;
sparsity_params.eta=0.1;
sparsity_params.sparsity=1;

[gamma_path mu_path sigma_path total_time soft_assignments bg_rate]= ...
    EM_fullmodel_v3(mpp(1:n_trial), ...
    event_rates,...
    evoked_cell_stimulated,expected_all, ...
    n_cell_local, gamma_old, mu_old, sigma_old, ...
    convergence_epsilon,f_background, mean_background, sigma_background, ...
    maxit,t_vect,use_size,background_update,sparsity_params);

gamma_ini= gamma_path(:,end);
sigma_ini = sigma_path(:,end);
mu_ini = mu_path(:,end);

%% Visualize the initial fits
figure(1)

    gamma_initial =zeros(n_cell_local,1);
    gamma_initial(stimulated_cells)=gamma_ini;
         plot(local_gamma(stimulated_cells)+normrnd(0,0.01,[sum(stimulated_cells) 1]),gamma_initial(stimulated_cells),'.','MarkerSize',20)
     xlim([-0.1,1.1]);
ylim([-0.1,1.1]);

line([0 1],[0 1]);
xlabel('True gamma')
ylabel('Est. gamma')

% save('initial_fits_Set2.mat','gamma_initial');
