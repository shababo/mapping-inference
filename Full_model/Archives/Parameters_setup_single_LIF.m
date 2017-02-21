%% Parameter setups for the LIF
load('../Optimal_design/Environments/current_template.mat'); %Contains the vector norm_average_current

I_eg_100mV=abs(norm_average_current(1:20:20*75).*100);
I_eg_50mV=abs(norm_average_current(1:20:20*75).*50);
I_eg_25mV=abs(norm_average_current(1:20:20*75).*25);
I_e=[repmat(I_eg_100mV',1,5) repmat(I_eg_50mV',1,5) repmat(I_eg_25mV',1,5)];
%% New parameters 
stoc_mu=0;stoc_sigma=0.3;
g=0.1; %membrane time constant [ms]
k_offset = 0.004; % From the simulation

load('../Optimal_design//Environments/fits_old.mat'); %
%% generate cell features conditioned on location

% layer based priors on featues
cell_feature_priors.connection_prob = 0.3; % bernoulli
cell_feature_priors.connection_strength_mean = 4;
cell_feature_priors.connection_strength_stddev = connection_strength_stddev;

% Using values in the simulations
cell_feature_priors.Vthre_mean = -25;
cell_feature_priors.Vthre_std = Vthre_std;

% Need to set the resting potential to be lower than the Vthre_mean (for
% each neuron)
cell_feature_priors.Vreset_mean = -80;
cell_feature_priors.Vreset_std = Vreset_std;
% Difference between resting potential and the firing threshold:
cell_features_priors.dE_L= 5;

all_connected = rand(n_cell,1) < cell_feature_priors.connection_prob;
    
all_amplitudes = normrnd(cell_feature_priors.connection_strength_mean,...
        cell_feature_priors.connection_strength_stddev,...
        [n_cell 1]);
all_V_th = normrnd(cell_feature_priors.Vthre_mean,...
        cell_feature_priors.Vthre_std,...
        [n_cell 1]);
all_V_reset = normrnd(cell_feature_priors.Vreset_mean,...
        cell_feature_priors.Vreset_std,...
        [n_cell 1]);
    
all_E_L = all_V_th - cell_features_priors.dE_L;


%% voltage-clamp parameters (daq stuff, bg psc parameters)


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
evoked_params.failure_prob = 0.2; % 
