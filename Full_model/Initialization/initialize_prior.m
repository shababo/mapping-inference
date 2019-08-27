%% Create the prior distributions base on the processed pilot data
function  [prior_info]=initialize_prior(params)
%% Initialize the output:
prior_info=struct;
prior_info.prior_parameters=struct;
prior_info.background_rate=params.background_rate;

%% Load and process the prior distributions from the pilot data 
prior_params=struct;
prior_params.gain.dist='logit-normal';
prior_params.gain.type='individual';
prior_params.gain.mean=params.initial_mean;
prior_params.gain.log_sigma=log(params.initial_var);
prior_params.gain.bounds.up=params.bounds.gain(2);
prior_params.gain.bounds.low=params.bounds.gain(1);

if isfield(params,'gain_pilot_path')
    load(params.gain_pilot_path)
    fields_list=fieldnames(gain_pilot);
    % Define the centers as the point on the grid that has the minimal squared
    % distance to the origin
    gain=[];
    for i_fld = 1:length(fields_list)
        fld=fields_list{i_fld};
        %size(rmmissing(gain_pilot.(fld).gain))
        gain=[gain rmmissing(gain_pilot.(fld).gain)];
    end
    prior_params.gain.bounds.up=max(params.bounds.gain(2),max(gain)*1.1); 
    prior_params.gain.bounds.low=min(params.bounds.gain(1),min(gain)*(1-0.1));
    % logit transformation for given the lower and upper bound:
    logit_gain=log((gain-prior_params.gain.bounds.low)./(prior_params.gain.bounds.up-gain));
    prior_params.gain.mean=mean(logit_gain);
    prior_params.gain.log_sigma=log(var(logit_gain)); % This is actually the variance!
end

% Pair-patched data needed for these: 
prior_params.delay_mean.dist='logit-normal';
prior_params.delay_mean.type='individual';
prior_params.delay_mean.mean=params.initial_mean;
prior_params.delay_mean.log_sigma=log(params.initial_var);
prior_params.delay_mean.bounds.up=params.bounds.delay_mean(2);
prior_params.delay_mean.bounds.low=params.bounds.delay_mean(1);

prior_params.delay_var.dist='logit-normal';
prior_params.delay_var.type='individual';
prior_params.delay_var.mean=params.initial_mean;
prior_params.delay_var.log_sigma=log(params.initial_var);
prior_params.delay_var.bounds.up=params.bounds.delay_var(2);
prior_params.delay_var.bounds.low=params.bounds.delay_var(1);

% Not available from pilot: 
prior_params.PR.dist='logit-normal';
prior_params.PR.type='individual';
prior_params.PR.mean=params.initial_mean;
prior_params.PR.log_sigma=log(params.initial_var);
prior_params.PR.bounds.up=params.bounds.PR(2);
prior_params.PR.bounds.low=params.bounds.PR(1);
%         variational_params.PR.bounds.prob_logit=;
prior_params.background=struct;
prior_params.background.dist='logit-normal';
prior_params.background.mean=0;
prior_params.background.log_sigma=1;
prior_params.background.bounds=struct;
prior_params.background.bounds.up=2e-4;
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



prior_info.prior_parameters=prior_params;

%% Handle the GPs 
load(params.gp_pilot_path) % An object called gp_pilot
axis_list = fieldnames(gp_pilot);
% Add an additional variance term in the GP variance:
for i_ax = 1:length(axis_list)
    ax=axis_list{i_ax};
    if isfield(ax, params.GP_params.GP_var_additional)
        gp_pilot.(ax).var_params.values= gp_pilot.(ax).var_params.values+params.GP_params.GP_var_additional.(ax);
    end
end
prior_info.GP_params=gp_pilot; % Same infor regarding prior distributions in the subfield prior_parameters
prior_info.GP_params.GP_minimal_variance = params.GP_params.GP_minimal_variance;
prior_info.GP_params.GP_added_variance=params.GP_params.GP_added_variance;
prior_info.GP_params.GP_boundary=params.GP_params.GP_boundary;
prior_info.GP_params.type= params.GP_params.type;
prior_info.GP_params.min_dist= params.GP_params.min_dist;

%% Real the relationships between the current and mean spike time and variance spike time
load(params.spike_curve_path);
prior_info.induced_intensity=induced_intensity;
prior_info.induced_intensity.minimum_stim_threshold = ...
    prior_info.induced_intensity.current(max(find(prior_info.induced_intensity.prob<0.2))+1);
prior_info.induced_intensity.stim_threshold=params.induced_intensity.stim_threshold;

prior_info.induced_intensity.spike_time_max=params.induced_intensity.spike_time_max;
prior_info.induced_intensity.event_time_max=params.induced_intensity.event_time_max;
prior_info.induced_intensity.time_factor=params.time_factor;

%% Use pilot data to set the priors for these parameters:
% Use prior info, if path to the summary data is available
% Use set parameters, if path is not available


%% Examine the shifts (hardly any in the gp data)
% There is hardly any shift in these data...

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