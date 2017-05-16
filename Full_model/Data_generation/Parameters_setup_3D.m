%% Parameter setups for the LIF-GLM model
%% Parameters for the LIF model (use the new fits)
% Stochastic components of voltages 
stoc_mu=0;stoc_sigma=0.5; % so that we can allow for a coarse grid..

%% voltage-clamp parameters (daq stuff, bg psc parameters)
bg_params.mean = 0;
bg_params.sigma = 1.5;
bg_params.firing_rate = 10; %spike/sec 
%% Generate data from circuit mapping model with multi-cell stimuli
%% build layers
% set priors on layer boundaries (from Lefort et al 2009)

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
% exc_neurons_per_layer_mean = [0 546 1145 1656 454 641 1288]/3;
% exc_neurons_per_layer_sd = [0 120 323 203 112 122 205]/3;
exc_neurons_per_layer_mean = [0 546 1145 1656 454 641 1288]/9;
exc_neurons_per_layer_sd = [0 120 323 203 112 122 205]/9;
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
% Note: We should set the prioirs for Vthres and Vreset according to our
% experiments
% The synaptic success rates and amplitude distributions 
cell_feature_priors.connection_prob = [0 .095 .057 .116 .191 .017 .006]; % bernoulli
cell_feature_priors.connection_strength_mean = exp([0 .8 .6 .8 2.0 .4 .1]); % log-normal
%cell_feature_priors.connection_strength_mean = exp([0 .8 .6 .8 1.2 .4 .1]+1.2); % log-normal
cell_feature_priors.connection_strength_stddev = 0.5*ones(num_layers,1); % log-normal
% Spiking thresholds

% draw features for each neuron
neuron_features = struct();
for i = 1:num_layers
    num_neurons_layer = size(neuron_locations{i},1);
    % Draw the connectivity 
    neuron_features(i).connected = ...
        rand(num_neurons_layer,1) < cell_feature_priors.connection_prob(i);
    % mean and variance for the amplitude distribution
    neuron_features(i).amplitude = abs(normrnd(cell_feature_priors.connection_strength_mean(i),...
        cell_feature_priors.connection_strength_stddev(i),...
        [num_neurons_layer 1]));
    neuron_features(i).amplitude = neuron_features(i).amplitude .* neuron_features(i).connected;
    
    if trials_specific_variance == 1
    neuron_features(i).sigma_across = sqrt(0.5)*ones(num_neurons_layer,1);
    neuron_features(i).sigma_within = sqrt(0.5)*ones(num_neurons_layer,1);
    else
    neuron_features(i).sigma_across = 0*ones(num_neurons_layer,1);
    neuron_features(i).sigma_within = 1*ones(num_neurons_layer,1);
    end
    % Successful probability: positively propotional to the mean amplitude
    neuron_features(i).success_prob = neuron_features(i).connected.*...
        exp(neuron_features(i).amplitude/2)./(exp(neuron_features(i).amplitude/2)+1);
    % Draw the spiking thresholds for neurons in each layer
    neuron_features(i).V_th = normrnd(cell_feature_priors.Vthre_mean(i),...
        cell_feature_priors.Vthre_std(i),...
        [num_neurons_layer 1]);
    % Draw the reseting threshold for neurons in each layer
    neuron_features(i).V_reset = normrnd(cell_feature_priors.Vreset_mean(i),...
        cell_feature_priors.Vreset_std(i),...
        [num_neurons_layer 1]);
    neuron_features(i).V_th_layer = cell_feature_priors.Vthre_mean(i)*ones(num_neurons_layer,1);
    neuron_features(i).V_reset_layer = cell_feature_priors.Vreset_mean(i)*ones(num_neurons_layer,1);
    
    neuron_features(i).shape_gain = randsample(1:num_types_cell,num_neurons_layer,true);
    
end
%% condense all cells into single arrays (for data generation)
all_locations = [];
all_amplitudes = [];
all_sigma_within = [];
all_sigma_across = [];
all_V_th = [];
all_V_reset = [];
all_V_th_layer = [];
all_V_reset_layer = [];
all_gamma = [];
all_layer=[];
all_shape_gain = [];
for i = 1:num_layers
    all_locations = [all_locations; neuron_locations{i}];
    all_layer = [all_layer; i*ones(size(neuron_locations{i},1),1)];
    all_amplitudes = [all_amplitudes; neuron_features(i).amplitude];
    all_sigma_within = [all_sigma_within; neuron_features(i).sigma_within];
    all_sigma_across = [all_sigma_across; neuron_features(i).sigma_across];
    
    all_V_th = [all_V_th; neuron_features(i).V_th];
    all_V_reset = [all_V_reset; neuron_features(i).V_reset];
    all_V_th_layer = [all_V_th_layer; neuron_features(i).V_th_layer];
    all_V_reset_layer = [all_V_reset_layer; neuron_features(i).V_reset_layer];
    all_gamma = [all_gamma; neuron_features(i).success_prob];
    all_shape_gain = [all_shape_gain; neuron_features(i).shape_gain'];
end
all_connected = all_amplitudes>0;
n_cell = length(all_amplitudes); % num_neurons


%% select a postsyanptic cell
cell_layer = 3; % layer 3
num_cell_layer_neurons = size(neuron_locations{cell_layer},1);
postsyn_position = zeros(1,3);
while postsyn_position(1) < 100 || postsyn_position(1) > 500
    postsyn_position = neuron_locations{cell_layer}(randi(num_cell_layer_neurons),:);
end


%%
region_width = 300;
region_height = 300;

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

% Identify neurons within this region:
neuron_in_region = zeros(size(all_neuron_locations,1),1); 
for i = 1:size(all_neuron_locations,1)
    if all_neuron_locations(i,1) > x_low & all_neuron_locations(i,1) < x_upp
        if all_neuron_locations(i,2) > y_low & all_neuron_locations(i,2) < y_upp
            neuron_in_region(i)=1;
        end
    end
end

%% Save info of the local neurons (for inference)
local_locations = [];
local_amplitudes = [];
local_V_th = [];
local_sigma_within= [];
local_sigma_across= [];
local_V_reset = [];
local_gamma = [];
local_index = [];
local_V_th_layer = [];
local_V_reset_layer = [];
local_layer=[];
local_shape_gain=[];
for i = 1:n_cell
     if neuron_in_region(i) > 0
   local_index = [local_index; i];
   local_layer = [local_layer; all_layer(i)];
   
    local_locations = [local_locations; all_locations(i,:)];
    local_amplitudes = [local_amplitudes; all_amplitudes(i)];
    local_V_th = [local_V_th; all_V_th(i)];
    local_V_reset = [local_V_reset; all_V_reset(i)];
     local_V_th_layer = [local_V_th_layer; all_V_th_layer(i)];
    local_V_reset_layer = [local_V_reset_layer; all_V_reset_layer(i)];
    local_sigma_within = [local_sigma_within; all_sigma_within(i)];
    local_sigma_across = [local_sigma_across; all_sigma_across(i)];
    local_gamma = [local_gamma; all_gamma(i)];
    local_shape_gain = [local_shape_gain; all_shape_gain(i)];
     end
end
local_connected = local_amplitudes>0;
Z = local_locations;
n_cell_local = size(Z,1); % Number of neurons in the region


%% Stimulus locations:
% For now, use only the nuclei

% Define a num_dense by num_dense by num_dense grid
% buffer=50;
% x_dense = (0:(num_dense-1))*(2*buffer+max(Z(:,1))-min(Z(:,1)))/(num_dense-1) + min(Z(:,1))-buffer;
% y_dense = (0:(num_dense-1))*(2*buffer+max(Z(:,2))-min(Z(:,2)))/(num_dense-1) + min(Z(:,2))-buffer;
% z_dense = (0:(num_dense-1))*(2*buffer+max(Z(:,3))-min(Z(:,3)))/(num_dense-1) + min(Z(:,3))-buffer;
% Z_dense = zeros(num_dense^2,3);
% for i = 1:num_dense
%     for l = 1:num_dense
%         for k = 1:num_dense
%         Z_dense((i-1)*num_dense + l,:) = [x_dense(i) y_dense(l) z_dense(k)];
%         end
%     end
% end
% % Merge the grid with the centers of neurons 
% Z_dense = [Z_dense; Z];

Z_dense = Z;

%% Assign shapes to local neurons using the templates
% Each shape is saved as a n_grid_shape^3 array describing the 200^3 um
% space centerd around the neuron

% square_radius = 20;
% shape_grid_n = 21;
% shape_grid_gap = square_radius/(shape_grid_n-1);
% shape_grid_x = linspace(-square_radius,square_radius,shape_grid_n);
% shape_grid_y = linspace(-square_radius,square_radius,shape_grid_n);
% shape_grid_z = linspace(-square_radius,square_radius,shape_grid_n);
% 
% all_shape_grid = zeros(n_cell,shape_grid_n,shape_grid_n,shape_grid_n);
% all_shape_mat = cell(n_cell,1);
% for i = 1:n_cell
%     radius_temp=radius{all_layer(i)};
%     i_ind = 0;
%     all_shape_mat{i} = zeros(0,3);
%     for i_x = 1:shape_grid_n
%         for i_y = 1:shape_grid_n
%             for i_z = 1:shape_grid_n
%                 distance_shape = sqrt(shape_grid_x(i_x)^2+shape_grid_y(i_y)^2+shape_grid_z(i_z)^2);
%                 if distance_shape < radius_temp
%                    all_shape_grid(i,i_x,i_y,i_z)=1;
%                    i_ind=i_ind+1;
%                    all_shape_mat{i}(i_ind,:) = all_locations(i,:) + [shape_grid_x(i_x) shape_grid_y(i_y) shape_grid_z(i_z)]; 
%                 end
%             end
%         end
%     end 
% end

