%% Generate data from circuit mapping model with multi-source stimuli
% 
% set RNG seed
rng(12242,'twister');

% Parameters on the grid and experimental design 
num_sources = 4; 

num_all = 20; 
num_grids = 31; 
grid_spacing = 15;

num_repeats = 1; % number of replicates

design=2;
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


%% Stimulating multiple spots in each trial


% these parameters govern the time delay, as a function of the
% point-spread-function stimuli for a particular trial
% in seconds
d_mean0 = .000;
d_sigma0 = .002;
d_mean_coef = .005;
d_sigma_coef = .050;

% We define the Z array as the locations of light stimuli
% This array contains the location of the stimuli inputs
% The dimension of array is N
Z = zeros(num_grids*num_grids,3);
grid_locations = zeros(num_grids*num_grids,2);
grid_index = zeros(num_grids,num_grids); 
for i = 1:num_grids
    for j = 1:num_grids
        grid_locations((i-1)*num_grids + j,:) = [i j];
        grid_index(i,j) = (i-1)*num_grids + j;
        Z((i-1)*num_grids + j,1) = (i-1)*grid_spacing - grid_spacing*floor(num_grids/2) + postsyn_position(1);
        Z((i-1)*num_grids + j,2) = (j-1)*grid_spacing - grid_spacing*floor(num_grids/2) + postsyn_position(2);
        Z((i-1)*num_grids + j,3) = postsyn_position(3); 
    end
end

%% Define the combination of trials 


% Option 1: almost pure randomization
% We define the combination by randomly permuting the sequence 1:num_grids^2,
% and then let num_sources consequent spots to be stimulated in each trial.
% For the last trial, we enroll elements at the beginning of the sequence
% if it is not full. 
if design == 1
    num_combinations = ceil(num_grids*num_grids/num_sources)*num_all;
    trial_locations_on_grid = zeros(num_combinations, num_sources);
    for m = 1:num_all
        perm_sequence = randperm(num_grids*num_grids);
        for	i = 1:ceil(num_grids*num_grids/num_sources)
            if i < ceil(num_grids*num_grids/num_sources)
                trial_locations_on_grid(i + (m-1)*num_combinations/num_all,:) = perm_sequence( (i-1)*num_sources + (1:num_sources));
            else
                num_missing = num_grids*num_grids-(i-1)*num_sources;
                trial_locations_on_grid(i + (m-1)*num_combinations/num_all,:) = [perm_sequence( ((i-1)*num_sources):end) perm_sequence(1:2)];
            end
        end
    end
else
    
    % Option 2: stimulate distant sites
    % To avoid singularity in design (the stimulated sets have rankds smaller
    % than the number of sites), we need an extra layer of randomness by alternating
    % between evenly spaced rows and columns
    gap=floor(num_grids/num_sources);
    num_combinations = num_grids*(num_grids- num_sources*gap+gap)*num_all;
    num_runs = num_combinations*num_repeats;
    trial_locations_on_grid = zeros(num_combinations, num_sources);
    for m = 1:num_all
        % Alternating between rows and columns!
        if mod(m,2)==0
            perm_index = grid_index;
            for i = 1:num_grids  % permuting the columns
                perm_index(:,i) =  perm_index(randperm(num_grids),i);
            end
            for i = 1:(num_grids- num_sources*gap+gap)
                col_index = (1:num_sources)*gap+i-gap;
                if max(col_index)>num_grids
                    col_index = mod(col_index,num_grids);
                end
                trial_locations_on_grid( (1:num_grids)+ (i-1)*num_grids + (m-1)*num_combinations/num_all,:) = perm_index(:,col_index);
            end
        else
            perm_index = grid_index;
            for i = 1:num_grids  % permuting the columns
                perm_index(i,:) =  perm_index(i,randperm(num_grids));
            end
            for i = 1:(num_grids- num_sources*gap+gap)
                row_index = (1:num_sources)*gap+i-gap;
                if max(row_index)>num_grids
                    row_index = mod(row_index,num_grids);
                end
                trial_locations_on_grid((1:num_grids)+ (i-1)*num_grids + (m-1)*num_combinations/num_all,:) = perm_index(row_index,:)';
            end
        end
    end
    
end
% covariates = zeros(size(trial_locations_on_grid,1), num_grids^2);
% for i = 1:num_combinations
% 	covariates(i, trial_locations_on_grid(i,:)) = 1;    
%  end
%  for i = 1:ceil(num_grids/num_sources)
%     col_index = (1:num_sources)*gap+i-gap
%  end
% rank(covariates)
% size(covariates)

%% Calculate the light-induced probability 
trial_locations_on_grid = repmat(trial_locations_on_grid,num_repeats,1);
N= size(trial_locations_on_grid,1);
pi_k = zeros(N,size(all_locations,1));

B = diag(all_locations*inv(A)*all_locations')*ones(1,num_sources);

for n = 1:N
   this_trial = trial_locations_on_grid(n,:);
   this_trial_locations = Z(this_trial,:);
   
   % for r = 1:R
	   % Pi(n,:,r) = exp(-sum(((p - p(:,stimulus(r)*ones(K,1)))'/A).*(p - p(:,stimulus(r)*ones(K,1)))',2)/2);
   % end
   % Calculate the differences between pairs of locations 
   % pi_kr = exp(-0.5*squareform(pdist([this_trial_locations; all_locations],'mahalanobis',A)).^2);
   
   diff_mat=B+(ones(size(all_locations,1),1)*diag(this_trial_locations*inv(A)*this_trial_locations')')-2*all_locations*inv(A)*this_trial_locations';

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





