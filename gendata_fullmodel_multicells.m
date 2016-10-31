%% Generate data from circuit mapping model with multi-cell stimuli
% NOTE: the locations of neurons are assumed to be known
% set RNG seed
rng(12242,'twister');

% Parameters on the grid and experimental design 
num_sources = 4; 
N = 6000;

% The map will be a square centered at the patched neuron
region_width = 500;
region_height = 500;

%% Gen neurons and their types/locations/features

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

% plot layers
% figure(1234)
% plot([0 1],-bsxfun(@times,[ones(num_layers + 1,2)],[layer_boundaries]'))

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

% layer based priors on featues
cell_feature_priors.connection_prob = [0 .095 .057 .116 .191 .017 .006]; % bernoulli
cell_feature_priors.connection_strength_mean = [0 .8 .6 .8 2.0 .4 .1]; % log-normal
cell_feature_priors.connection_strength_stddev = ones(num_layers,1); % log-normal
cell_feature_priors.connection_tau_rise_mean = [0 2.8 1.86 1.77 2.37 5.41 1.1]/1/1000; % gaussian
cell_feature_priors.connection_tau_rise_std = .0005*ones(num_layers,1)/2.5; % gaussian
cell_feature_priors.connection_tau_fall_mean = [0 73.2 37.2 61.7 74.4 36.8 27.2]/20/1000; % gaussian
cell_feature_priors.connection_tau_fall_std = .003*ones(num_layers,1)/10; % gaussian
cell_feature_priors.rheobase_mean = [0 126 132 56 68 98 76]/126; % gaussian
cell_feature_priors.rheobase_std = 5 * ones(num_layers,1); % gaussian

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
                                       
    neuron_features(i).tau_rise = normrnd(cell_feature_priors.connection_tau_rise_mean(i),...
                                           cell_feature_priors.connection_tau_rise_std(i),...
                                           [num_neurons_layer 1]);
    neuron_features(i).tau_rise(neuron_features(i).tau_rise < 0) = .001;                                   
                                       
    neuron_features(i).tau_fall = normrnd(cell_feature_priors.connection_tau_fall_mean(i),...
                                           cell_feature_priors.connection_tau_fall_std(i),...
                                           [num_neurons_layer 1]);    
                                       
	neuron_features(i).rheobase = normrnd(cell_feature_priors.rheobase_mean(i),...
                                           cell_feature_priors.rheobase_std(i),...
                                           [num_neurons_layer 1]); 
end

%% condense all cells into single arrays
all_locations = [];
all_amplitudes = [];
all_tau_rise = [];
all_tau_fall = [];
all_rheobase = [];

for i = 1:num_layers
    all_locations = [all_locations; neuron_locations{i}];
    all_amplitudes = [all_amplitudes; neuron_features(i).amplitude];
    all_tau_rise = [all_tau_rise; neuron_features(i).tau_rise];
    all_tau_fall = [all_tau_fall; neuron_features(i).tau_fall];
    all_rheobase = [all_rheobase; neuron_features(i).rheobase];
end

K = length(all_amplitudes); % num_neurons



%% voltage-clamp parameters (daq stuff, bg psc parameters)
tau_r_bounds = [1 20]/20000;
tau_f_bounds = [50 200]/20000;

data_params.T = 2000; %bins - start not too long
data_params.dt = 1/20000; %
data_params.baseline = 0;
data_params.sigmasq = 3.5;
data_params.phi = [1, .80, -.12]; %this determines what the AR noise looks like.


bg_params.tau_r_bounds = tau_r_bounds*20000;
bg_params.tau_f_bounds = tau_f_bounds*20000;
bg_params.a_min = .5;
bg_params.a_max = 15;
bg_params.firing_rate = 20; %spike/sec 

%% stim paramters

% covariance of point spread function
A = diag([200, 200, 750]);

evoked_params.stim_tau_rise = .0015*20000; % values for chr2 from lin et al 2009 (biophysics)
evoked_params.stim_tau_fall = .013*20000;
evoked_params.stim_amp = 0;
evoked_params.stim_start = .005*20000;
evoked_params.stim_duration = .005*20000;

% effect on postsynaptic cell
evoked_params.stim_tau_rise = .0015*20000; % values for chr2 from lin et al 2009 (biophysics)
evoked_params.stim_tau_fall = .013*20000;
evoked_params.stim_amp = 0;

%% select a postsyanptic cell
cell_layer = 5; % 5A
num_cell_layer_neurons = size(neuron_locations{cell_layer},1);

postsyn_position = zeros(1,3);
while postsyn_position(1) < 100 || postsyn_position(1) > 500
    postsyn_position = neuron_locations{cell_layer}(randi(num_cell_layer_neurons),:);
end


%% Simulating Experiments:

%% Find the neurons that locate in the region:
% 
% all_neuron_locations record locations of the neurons and their layer info
all_neuron_locations =zeros(0,4);
for i = 1:num_layers
    all_neuron_locations = [all_neuron_locations; [neuron_locations{i}(:,1:3) i*ones(size(neuron_locations{i},1),1) ]];
end

epsilon = 0.01;
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

Z = zeros(sum(neuron_in_region),3);
count = 1;
for i = 1:size(all_neuron_locations,1)
    if neuron_in_region(i) > 0
        Z(count,:) = all_neuron_locations(i,1:3);
        count = count + 1;
    end 
end


%% Define the combination of trials 
% We want to stimulate neurons that are far apart from each other

% To do this, we split the cells by the quantiles of their x and y
% coordinates, then alternating between x and y.

% Find the quantiles:
nquantile = num_sources*2;

[~, ~, x_ranks] = unique(Z(:,1));
x_freq = x_ranks / size(Z,1);
[~, ~, y_ranks] = unique(Z(:,2));
y_freq = y_ranks / size(Z,1);

x_index = ceil(x_freq / (1/nquantile));
y_index = ceil(y_freq / (1/nquantile));

x_group = cell(nquantile,1);
y_group = cell(nquantile,1);
for i = 1:size(Z,1)
    x_group{x_index(i)} = [x_group{x_index(i)} i];
    y_group{y_index(i)} = [y_group{y_index(i)} i];
end
max_x_group = 0;
max_y_group = 0;
for i = 1:nquantile
    max_x_group = max(max_x_group, size(x_group{i},2));
    max_y_group = max(max_y_group, size(y_group{i},2));
end

joint_n = 2*(max_x_group + max_y_group); 
M = ceil(N/joint_n);

trial_locations_on_grid = zeros(N, num_sources);
count = 0;
for m = 1:M
    % Alternating between rows and columns!
    
    perm_index = ones(max_x_group,nquantile);
    for i = 1:nquantile  % permuting the columns
        n_this_group = size(x_group{i},2);
        perm_index(1:n_this_group,i) = x_group{i}(randperm(size(x_group{i},2)));
    end
    trials_m = [perm_index(:, 2*(1:num_sources)-1); perm_index(:, 2*(1:num_sources))];
    
    perm_index = ones(max_y_group,nquantile);
    for i = 1:nquantile  % permuting the columns
        n_this_group = size(y_group{i},2);
        perm_index(1:n_this_group,i) =  y_group{i}(randperm(size(y_group{i},2)));
    end
    trials_m = [trials_m; [perm_index(:, 2*(1:num_sources)-1); perm_index(:, 2*(1:num_sources))]];
    
    
    if N < (m*joint_n) 
        trial_locations_on_grid( ((m-1)*joint_n+1 ):N,:) = trials_m(1: (N-((m-1)*joint_n )),:);
    else
        trial_locations_on_grid( (m-1)*joint_n + (1:joint_n),:) = trials_m;
    end 
end

% covariates = zeros(size(trial_locations_on_grid,1), size(Z,1));
% for i = 1:N
% 	covariates(i, trial_locations_on_grid(i,:)) = 1;    
% end
% rank(covariates)
% size(covariates)

%% Calculate the light-induced probability 
% these parameters govern the time delay, as a function of the
% point-spread-function stimuli for a particular trial
% in seconds
d_mean0 = .000;
d_sigma0 = .002;
d_mean_coef = .005;
d_sigma_coef = .050;

pi_k = zeros(N,size(all_locations,1));
B = diag(all_locations(:,1:3)*inv(A)*all_locations(:,1:3)')*ones(1,num_sources);

for n = 1:N
   this_trial = trial_locations_on_grid(n,:);
   this_trial_locations = Z(this_trial,:);
   
   % for r = 1:R
	   % Pi(n,:,r) = exp(-sum(((p - p(:,stimulus(r)*ones(K,1)))'/A).*(p - p(:,stimulus(r)*ones(K,1)))',2)/2);
   % end
   % Calculate the differences between pairs of locations 
   % pi_kr = exp(-0.5*squareform(pdist([this_trial_locations; all_locations],'mahalanobis',A)).^2);
   
   diff_mat=B+(ones(size(all_locations,1),1)*diag(this_trial_locations*inv(A)*this_trial_locations')')-2*all_locations(:,1:3)*inv(A)*this_trial_locations';

   pi_kr = exp(-0.5*(diff_mat));
   pi_k(n,:) = min(.95,sum(pi_kr,2));
end

pi_k_spike = pi_k;
pi_k_spike(pi_k_spike > .65) = 1; % what does this mean?

% % firing delay means and variances
% d_mean_nk = d_mean0 + (1 - pi_nk)*d_mean_coef;
% d_sigma_nk = d_sigma0 + (1 - pi_nk)*d_sigma_coef;
% 
% % sample "ground truth" firing delay
% D = normrnd(d_mean_nk,d_sigma_nk)/data_params.dt + evoked_params.stim_start;
% D(D < evoked_params.stim_start + .002) = evoked_params.stim_start + .002;

% firing delay means and variances
d_mean_nk = d_mean0 + (1.5 - pi_k)*d_mean_coef;
d_sigma_nk = d_sigma0 + (1 - pi_k)*d_sigma_coef;

% sample "ground truth" firing delay
D = exprnd(d_mean_nk/data_params.dt) + evoked_params.stim_start;
% D(D < evoked_params.stim_start + .002) = evoked_params.stim_start + .002;

% sample "ground truth" stimulations

X = rand(N,K) < pi_k_spike; %.2 
X(D > 2000) = 0;

%% Generate a response, given D, Pi, X, w

Y = zeros(N,data_params.T);

for n = 1:N
    firing_neurons = X(n,:) & all_amplitudes' > 0;
    
    if any(all_amplitudes(firing_neurons) > 0)
        
        evoked_params.times = D(n,firing_neurons);
        evoked_params.a = all_amplitudes(firing_neurons);
        evoked_params.tau_r = all_tau_rise(firing_neurons)/data_params.dt;
        evoked_params.tau_f = all_tau_fall(firing_neurons)/data_params.dt;
    else
        evoked_params.times = [];
        evoked_params.a = [];
        evoked_params.tau_r = [];
        evoked_params.tau_f = [];
    end
    
    [Y(n,:), mpp_n] = gen_trace(data_params,bg_params,evoked_params);
    if n == 1
        mpp = mpp_n;
    else
        mpp(n) = mpp_n;
    end
end





