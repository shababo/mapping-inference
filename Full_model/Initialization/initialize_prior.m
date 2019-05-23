%% Create the prior distributions base on the processed pilot data
function  [prior_info]=initialize_prior(gp_pilot_path,spike_curve_path)

%% Handle the GPs 
load(gp_pilot_path) % An object called gp_pilot
%% Examine the shifts (hardly any in the gp data)
% There is hardly any shift in these data 
% mean([gp_pilot.x.neurons(:).initial_shift])
% mean([gp_pilot.y.neurons(:).initial_shift])
% mean([gp_pilot.z.neurons(:).initial_shift])
% gp_pilot.xy.neurons(:).initial_shift

% prior_params.shift_x.dist='logit-normal';
% prior_params.shift_x.type='individual';
% prior_params.shift_x.mean=summary_results.x.shift_params.mean;
% prior_params.shift_x.log_sigma=log(summary_results.x.shift_params.var)/2;
% 
% prior_params.shift_y.dist='logit-normal';
% prior_params.shift_y.type='individual';
% prior_params.shift_y.mean=summary_results.y.shift_params.mean;
% prior_params.shift_y.log_sigma=log(summary_results.y.shift_params.var)/2;
% 
% prior_params.shift_z.dist='logit-normal';
% prior_params.shift_z.type='individual';
% prior_params.shift_z.mean=summary_results.z.shift_params.mean;
% prior_params.shift_z.log_sigma=log(summary_results.z.shift_params.var)/2;

%% Convolute the GP mean functions with an additional normal noise:

axis_list = fieldnames(pilot_data);
for i_ax = 1:length(axis_list)
    ax=axis_list{i_ax};
    neurons=pilot_data.(ax).neurons;
    n_cell = length(neurons);

end



%% Add an additional variance term in the GP variance:


%% 
load(spike_curve_path) % 


sum_path = [prior_path];
load(sum_path)

bounds=struct;
% bounds.PR=[0.01 1];
% bounds.gain=[0.005 0.06];
% bounds.delay_mu=[0 60];
% bounds.delay_sigma=[0.1 10];
bounds.PR=[0.01 1];
bounds.gain=[0.01 0.05];
bounds.delay_mu=[30 70];
bounds.delay_sigma=[5 10];

ini_mean=0;
ini_sigma=2;
prior_params=struct;

prior_params.GP_params=summary_results;
prior_params.GP_minimal_variance = 0.01;
%% Initialize the priors for all other parameters:
prior_params.boundary_params= [80 30 120];
prior_params.initial_boundary_params= [10 10 30];


prior_params.gain.dist='logit-normal';
prior_params.gain.type='individual';
prior_params.gain.mean=ini_mean;
prior_params.gain.log_sigma=log(ini_sigma);
prior_params.gain.bounds.up=bounds.gain(2);
prior_params.gain.bounds.low=bounds.gain(1);



prior_params.PR.dist='logit-normal';
prior_params.PR.type='individual';
prior_params.PR.mean=ini_mean;
prior_params.PR.log_sigma=log(ini_sigma);
prior_params.PR.bounds.up=bounds.PR(2);
prior_params.PR.bounds.low=bounds.PR(1);
%         variational_params.PR.bounds.prob_logit=;

prior_params.delay_mu.dist='logit-normal';
prior_params.delay_mu.type='individual';
prior_params.delay_mu.mean=ini_mean;
prior_params.delay_mu.log_sigma=log(ini_sigma);
prior_params.delay_mu.bounds.up=bounds.delay_mu(2);
prior_params.delay_mu.bounds.low=bounds.delay_mu(1);


prior_params.delay_sigma.dist='logit-normal';
prior_params.delay_sigma.type='individual';
prior_params.delay_sigma.mean=ini_mean;
prior_params.delay_sigma.log_sigma=log(ini_sigma);
prior_params.delay_sigma.bounds.up=bounds.delay_sigma(2);
prior_params.delay_sigma.bounds.low=bounds.delay_sigma(1);

prior_params.background=struct;
prior_params.background.dist='logit-normal';
prior_params.background.mean=0;
prior_params.background.log_sigma=1;
prior_params.background.bounds=struct;
prior_params.background.bounds.up=1e-2;
prior_params.background.bounds.low=1e-5;
prior_params.background.type='common';

prior_params.shapes=struct;
% prior_params.shapes.dist='logit-normal';
prior_params.shapes.dist='mvn';  % using GP approximation
prior_params.shapes.type='individual';
prior_params.shapes.locations=zeros(0,3);
prior_params.shapes.mean=zeros(0,1);
prior_params.shapes.log_sigma=zeros(0,1);
prior_params.shapes.bounds.up=zeros(0,1);
prior_params.shapes.bounds.low=zeros(0,1);
prior_params.shapes.prior_sigma=zeros(0,1);



prior_info.prior_parameters.GP_params.type='xy_square';
prior_info.induced_intensity=get_spike_curves(spike_curve_path);
prior_info.induced_intensity.minimum_stim_threshold = ...
    prior_info.induced_intensity.current(max(find(prior_info.induced_intensity.prob<0.2))+1);
prior_info.induced_intensity.stim_threshold=5;
prior_info.induced_intensity.spike_time_max=200;
prior_info.induced_intensity.event_time_max=300;
prior_info.background_rate=background_rate;
prior_info.prior_parameters.initial_boundary_params=[20 20 30];
if nuclei_design
    prior_info.prior_parameters.initial_boundary_params=[1 1 1];
end
prior_info.GP_added_variance=0.04;
