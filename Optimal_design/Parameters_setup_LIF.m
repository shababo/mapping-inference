%% Parameter setups for the LIF
load('./Environments/current_template.mat'); %Contains the vector norm_average_current

I_eg_100mV=abs(norm_average_current(1:20:20*75).*100);
I_eg_50mV=abs(norm_average_current(1:20:20*75).*50);
I_eg_25mV=abs(norm_average_current(1:20:20*75).*25);
I_e=[repmat(I_eg_100mV',1,5) repmat(I_eg_50mV',1,5) repmat(I_eg_25mV',1,5)];

%% New parameters 
%

%dt=1; %time step ms
% End of trial: 
%t_end=75;
% Stimulus onset time:
%t_StimStart=5;
%t_StimEnd=15;


% Stochastic components of voltages 
stoc_mu=0;stoc_sigma=0.3;


% Cell-specific values (need to draw new ones)
% V_th=-25; %spike threshold [mV]
% V_reset=-80; %value to reset voltage to after a spike [mV]

%E_L=-28; %resting membrane potential [mV]

g=0.1; %membrane time constant [ms]

k_offset = 0.004; % From the simulation

load('./Environments/fits_old.mat'); %
%mean(Vreset_g(2:end,1)) % -6.09
%std(Vreset_g(2:end,1)) % 0.786

%mean(Vreset_g(2:end,2)) % -143.6
%std(Vreset_g(2:end,2)) % 239.5

%gain_allLoc


%% Generate data from circuit mapping model with multi-cell stimuli
%% build layers
% set priors on layer boundaries (from Lefort et al 2009)
num_layers = 7;
layer_names = {'L1','L2','L3','L4','L5A','L5B','L6'};
layer_bottom_means = [0 128 269 418 588 708 890 1154];
layer_bottom_sds = [0 18 36 44 44 62 80 116]/10;
% draw layer boundaries
layer_boundaries = normrnd(layer_bottom_means,layer_bottom_sds);
while any(diff(layer_boundaries) < 20)
    layer_boundaries = normrnd(layer_bottom_means,layer_bottom_sds);
end
%% draw number of cells per layer
% set priors (from Lefort et al 2009)
exc_neurons_per_layer_mean = [0 546 1145 1656 454 641 1288]/3;
exc_neurons_per_layer_sd = [0 120 323 203 112 122 205]/3;
K_layers = ceil(normrnd(exc_neurons_per_layer_mean,exc_neurons_per_layer_sd));
while any(K_layers < 0)
    K_layers = ceil(normrnd(exc_neurons_per_layer_mean,exc_neurons_per_layer_sd));
end

% size of region containing neurons (or region we can stim)
barrel_width = 600;
slide_width = 100;

% how many neurons are excitatory (for now we will only consider a
% homogenous population of excitatory neurons - in the future we may
% consider many different cell types whose properties are different
pct_excitatory = 1.00;
neuron_locations = cell(num_layers,1);
for i = 1:num_layers
    neuron_locations{i} = sample_neuron_positions(K_layers(i), ...
        [0 barrel_width; layer_boundaries(i) layer_boundaries(i+1); 0 slide_width]);
end

%% generate cell features conditioned on location
% Note: 


% layer based priors on featues
cell_feature_priors.connection_prob = [0 .095 .057 .116 .191 .017 .006]; % bernoulli
cell_feature_priors.connection_strength_mean = [0 .8 .6 .8 2.0 .4 .1]; % log-normal
cell_feature_priors.connection_strength_stddev = ones(num_layers,1); % log-normal

% Using values in the simulations
cell_feature_priors.Vthre_mean = [0  -25  -25  -25  -25  -25  -25]; % gaussian
cell_feature_priors.Vthre_std = 0.2 * ones(num_layers,1); % gaussian

% Need to set the resting potential to be lower than the Vthre_mean (for
% each neuron)

cell_feature_priors.Vreset_mean = [0  -80  -80  -80  -80  -80  -80]; % gaussian
cell_feature_priors.Vreset_std = 16 * ones(num_layers,1); % gaussian

% Difference between resting potential and the firing threshold:
cell_features_priors.dE_L= [0  -3  -3  -3  -3  -3  -3]; % constant

% draw features for each neuron
neuron_features = struct();
for i = 1:num_layers
    
    num_neurons_layer = size(neuron_locations{i},1);
    
    neuron_features(i).connected = ...
        rand(num_neurons_layer,1) < cell_feature_priors.connection_prob(i);
    
    neuron_features(i).amplitude = lognrnd(cell_feature_priors.connection_strength_mean(i),...
        cell_feature_priors.connection_strength_stddev(i),...
        [num_neurons_layer 1]);
    neuron_features(i).amplitude = neuron_features(i).amplitude .* neuron_features(i).connected;
    
    %     neuron_features(i).tau_rise = normrnd(cell_feature_priors.connection_tau_rise_mean(i),...
    %                                            cell_feature_priors.connection_tau_rise_std(i),...
    %                                            [num_neurons_layer 1]);
    %     neuron_features(i).tau_rise(neuron_features(i).tau_rise < 0) = .001;
    %
    %     neuron_features(i).tau_fall = normrnd(cell_feature_priors.connection_tau_fall_mean(i),...
    %                                            cell_feature_priors.connection_tau_fall_std(i),...
    %                                            [num_neurons_layer 1]);
    
    neuron_features(i).rheobase = normrnd(cell_feature_priors.rheobase_mean(i),...
                                           cell_feature_priors.rheobase_std(i),...
                                           [num_neurons_layer 1]); 
                                       
    % New features needed in the LIF model
    neuron_features(i).V_th = normrnd(cell_feature_priors.Vthre_mean(i),...
        cell_feature_priors.Vthre_std(i),...
        [num_neurons_layer 1]);
    
            
    % New features needed in the LIF model
    neuron_features(i).E_L = neuron_features(i).V_th - cell_features_priors.dE_L(i);
    neuron_features(i).V_reset = normrnd(cell_feature_priors.Vreset_mean(i),...
        cell_feature_priors.Vreset_std(i),...
        [num_neurons_layer 1]);
    
end

%% condense all cells into single arrays
all_locations = [];
all_amplitudes = [];
%all_tau_rise = [];
%all_tau_fall = [];
all_rheobase = [];
all_V_th = [];
all_V_reset = [];
all_E_L = [];

for i = 1:num_layers
    all_locations = [all_locations; neuron_locations{i}];
    all_amplitudes = [all_amplitudes; neuron_features(i).amplitude];
    %all_tau_rise = [all_tau_rise; neuron_features(i).tau_rise];
    %all_tau_fall = [all_tau_fall; neuron_features(i).tau_fall];
    all_rheobase = [all_rheobase; neuron_features(i).rheobase];
    all_V_th = [all_V_th; neuron_features(i).V_th];
    all_V_reset = [all_V_reset; neuron_features(i).V_reset];
    all_E_L = [all_E_L; neuron_features(i).E_L];
end

K = length(all_amplitudes); % num_neurons
%% voltage-clamp parameters (daq stuff, bg psc parameters)
% tau_r_bounds = [1 20]/20000;
% tau_f_bounds = [50 200]/20000;


%dt=1; %time step ms
% End of trial: 
%t_end=75;
% Stimulus onset time:
%t_StimStart=5;
%t_StimEnd=15;

data_params.T = 75; % total time (ms)
data_params.dt = 1; % 1/20 ms

data_params.baseline = 0;
data_params.sigmasq = 3.5;
data_params.phi = [1, .80, -.12]; %this determines what the AR noise looks like.


% bg_params.tau_r_bounds = tau_r_bounds*20000;
% bg_params.tau_f_bounds = tau_f_bounds*20000;
bg_params.a_min = .5;
bg_params.a_max = 15;
bg_params.firing_rate = 10; %spike/sec 

%% stim paramters
evoked_params.stim_start = .005*1000/data_params.dt;
evoked_params.stim_duration = .010*1000/data_params.dt;

% effect on postsynaptic cell
% evoked_params.stim_tau_rise = .0015*20000; % values for chr2 from lin et al 2009 (biophysics)
% evoked_params.stim_tau_fall = .013*20000;
evoked_params.stim_amp = 0;

evoked_params.sigma_a = 0.2;
evoked_params.failure_prob = 0.1;
%% select a postsyanptic cell
cell_layer = 5; % 5A
num_cell_layer_neurons = size(neuron_locations{cell_layer},1);

postsyn_position = zeros(1,3);
while postsyn_position(1) < 100 || postsyn_position(1) > 500
    postsyn_position = neuron_locations{cell_layer}(randi(num_cell_layer_neurons),:);
end

%% Connectivity and amplitudes
% 
region_width = 500;
region_height = 500;

% all_neuron_locations record locations of the neurons and their layer info
all_neuron_locations =zeros(0,4);
for i = 1:num_layers
    all_neuron_locations = [all_neuron_locations; [neuron_locations{i}(:,1:3) i*ones(size(neuron_locations{i},1),1) ]];
end
epsilon = 0.01; % To avoid ties
x_low = postsyn_position(1) - region_width/2-epsilon;
x_upp = postsyn_position(1) + region_width/2;
y_low = postsyn_position(2) - region_height/2-epsilon;
y_upp = postsyn_position(2) + region_height/2;

% Identify neurons within this location:
neuron_in_region = zeros(size(all_neuron_locations,1),1); 
for i = 1:size(all_neuron_locations,1)
    if all_neuron_locations(i,1) > x_low & all_neuron_locations(i,1) < x_upp
        if all_neuron_locations(i,2) > y_low & all_neuron_locations(i,2) < y_upp
            neuron_in_region(i)=1;
        end
    end
end

Z = zeros(sum(neuron_in_region),3); % Stimulus locations
count = 1;
for i = 1:size(all_neuron_locations,1)
    if neuron_in_region(i) > 0
        Z(count,:) = [all_neuron_locations(i,1:2) postsyn_position(3)];
        count = count + 1;
    end 
end

K_z = size(Z,1); % Number of neurons in the region
 
local_neuron_amplitudes = zeros(K_z,1); % Amplitudes of neurons in this region

count = 1;
for i = 1:size(all_amplitudes,1)
    if neuron_in_region(i) > 0
        local_neuron_amplitudes(count) = all_amplitudes(i);
        count = count + 1;
    end 
end
