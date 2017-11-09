addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Load the data set 
load('./Environments/6_3_s2c3_mpp_and_stim_data.mat')
target_locations =target_locs;
single_trial_limit=max(find(isnan(target_inds(:,2))));
trials_locations= target_inds(1:single_trial_limit,:);
trials_powers = stim_pow(1:single_trial_limit);
load('./Environments/6_3_cell_locs_reduced_and_fixed.mat')
cell_locations = cell_locs;
n_trial = size(trials_locations,1);
n_cell=size(cell_locations,1);
%% Load the current template 
load('./Environments/chrome-template-3ms.mat');
downsamp=1;max_time=300;
power_level = unique(trials_powers);
num_power_level=length(power_level);
current_template=template(1:downsamp:max_time);
t_vect= 1:1:max_time;
%% Load the shape template 
load('./Environments/l23_template_cell.mat');
%l23_average_shape
temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;
shape_template=l23_average_shape;
%% Reduce the mpp to interested time period, and add the locations and powers 
min_time=0;
for i_trial = 1:n_trial
    mpp(i_trial).locations = trials_locations(i_trial,:);
    mpp(i_trial).power = trials_powers(i_trial,:);
    if mpp(i_trial).num_events >0
        range_idx = mpp(i_trial).times<max_time & mpp(i_trial).times>min_time ;
        mpp(i_trial).num_events = sum(range_idx);
        mpp(i_trial).times = mpp(i_trial).times(range_idx);
        mpp(i_trial).amp =mpp(i_trial).amp(range_idx);
    end
end
%%
% Estimate the background rate
v_reset_known=-4e3;
mpp_copy=mpp;
mpp=mpp(1:n_trial);
% Why are there a lot of events around 500?
trial_counts = 0;
event_counts = 0;
for i_trial = 1:n_trial
    %if stim_pow(i_trial)==50
    trial_counts = trial_counts+1;
    event_counts = event_counts +sum(mpp_copy(i_trial).times>20 & mpp_copy(i_trial).times<80);
    %end
end
background_rate = event_counts/trial_counts/60;
%% Paramters:
stim_threshold =10;g=0.02;v_th_known=15*ones(n_cell,1);
linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};

gain_grid=0.001*[8:30];gain_prior=normpdf(gain_grid,0.015,0.005);
gamma_grid= 0.1*[0:10]; 
gamma_prior=gamma_grid;gamma_prior(1)=0.7;
gamma_prior(2:end)= (1- gamma_prior(1))/(length(gamma_grid)-1);


n_gibbs_sample=100;

delay_params.type=2; %1: normal; 2: gamma
delay_params.mean=100; delay_params.std=30;
delay_params.delayed=true; delay_params.n_grid=200;
n_round_digit=0;
%% Initialization:
gain_initial= [];
gamma_initial = [];

% params.invlink = @invlink_sig;
% params.dlink = @derlink_sig;
% params.link = @link_sig;
% params.dinvlink = @derinvlink_sig;
% linkfunc = {params.link, params.dlink, params.invlink,params.dinvlink};
%% run the Gibbs sampler
[gain_samples, gamma_samples] = Gibbs_first_spike(mpp, ...
    target_locations, cell_locations,...
    current_template, shape_template,delay_params,linkfunc,...
    stim_threshold, g, background_rate,v_th_known,...
    gain_grid, gain_prior, gamma_grid,gamma_prior,...
    gamma_initial,gain_initial,n_gibbs_sample,...
     n_round_digit);
%%

% 1:
%delay_params.type=2; %1: normal; 2: gamma
%delay_params.mean=60; delay_params.std=15;
save('Gibbs_samples_6_3_July16.mat','gamma_samples','gain_samples');


% 2:
%delay_params.type=2; %1: normal; 2: gamma
%delay_params.mean=60; delay_params.std=30;
% save('./Results/Gibbs_samples_6_3_v2.mat','gamma_samples','gain_samples');


% 3:
%delay_params.type=2; %1: normal; 2: gamma
%delay_params.mean=35; delay_params.std=30;
% save('./Results/Gibbs_samples_6_3_v3.mat','gamma_samples','gain_samples');


% 4:
%delay_params.type=2; %1: normal; 2: gamma
%delay_params.mean=100; delay_params.std=30;
% save('./Results/Gibbs_samples_6_3_v4.mat','gamma_samples','gain_samples');





