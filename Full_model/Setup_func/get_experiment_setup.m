function experiment_setup = get_experiment_setup(varargin)



%experiment_setup.is_sim

% parameters

if ~isempty(varargin) && ~isempty(varargin{1})
    location_str = varargin{1};
else
    location_str = 'millennium_falcon';
end

% experiment; simulation; reproduction

experiment_setup.experiment_type='simulation'; % experiment; simulation; reproduction
    
switch location_str
    case 'szchen'
        experiment_setup.exp_root = 'C:/Users/Shizhe/Documents/Mapping_data/Data/';
        experiment_setup.analysis_root = 'C:/Users/Shizhe/Documents/Mapping_data/tmp/';   
        experiment_setup.experiment_type='simulation';
    case 'millennium_falcon'
        experiment_setup.exp_root = 'C:\data\Shababo';
        experiment_setup.analysis_root = '/media/shababo/data/'; % make sure to add ending slash
    case 'shababo'
        experiment_setup.exp_root = '/Users/shababo/projects/mapping/data/sim_tmp/';
        experiment_setup.analysis_root = '/Users/shababo/projects/mapping/data/sim_tmp/';
end
clock_array = clock;
experiment_setup.exp_id = [num2str(clock_array(2)) '_' num2str(clock_array(3)) ...
    '_' num2str(clock_array(4)) ...
    '_' num2str(clock_array(5))];
experiment_setup.exp.fullsavefile = fullfile(experiment_setup.exp_root,[experiment_setup.exp_id '_data.mat']);

% read this from argument
group_names=cell([4 1]);
group_names{1}='undefined';
group_names{2}='connected';
group_names{3}='disconnected';
group_names{4}='alive';
experiment_setup.group_names = group_names;
experiment_setup.default_group='undefined';

experiment_setup.terminator=@check_all_learned;

max_spots_per_trial = 0;
for i = 1:length(group_names)
    experiment_setup.groups.(group_names{i}) = eval(['get_' group_names{i}]);
    if isfield('design_func_params',experiment_setup.groups.(group_names{i}))
    max_spots_per_trial = max(max_spots_per_trial,...
        experiment_setup.groups.(group_names{i}).design_func_params.trials_params.spots_per_trial);
    end
end
experiment_setup.max_spots_per_trial = max_spots_per_trial;

% Get sim paramters
% sim params
experiment_setup.sim = get_simulation_setup;
 
switch location_str
    case 'szchen'
        experiment_setup.sim.compute_phase_masks=0;
        experiment_setup.sim.use_power_calib =0;
    case 'adesnik_lab'
        experiment_setup.sim.compute_phase_masks=1;
end





experiment_setup.prior_info=struct;

experiment_setup.prior_info.PR_prior = struct;
experiment_setup.prior_info.PR_prior.type='spiked_logit_normal'; % spike and logitnorm slab
experiment_setup.prior_info.PR_prior.pi_logit=-100;
experiment_setup.prior_info.PR_prior.alpha=-0.5;
experiment_setup.prior_info.PR_prior.beta=1;


experiment_setup.prior_info.gain_prior = struct;
experiment_setup.prior_info.gain_prior.type='spiked_logit_normal'; % spike and logitnorm slab
experiment_setup.prior_info.gain_prior.pi_logit=-Inf;
experiment_setup.prior_info.gain_prior.alpha=0;
experiment_setup.prior_info.gain_prior.beta=1;

%----------- Load the current template
load('chrome-template-3ms.mat');


experiment_setup.trials.downsamp = 1;
experiment_setup.trials.Fs = 20000;
experiment_setup.trials.min_time = 40;
experiment_setup.trials.max_time = 500;
experiment_setup.trials.max_time_sec = experiment_setup.trials.max_time/experiment_setup.trials.Fs;
current_template=template(1:experiment_setup.trials.downsamp:experiment_setup.trials.max_time);

experiment_setup.prior_info.current_template = current_template;
experiment_setup.prior_info.gain_model=[];
experiment_setup.prior_info.delay=struct;
experiment_setup.prior_info.delay.delayed=true;
experiment_setup.prior_info.delay.type='gamma';
experiment_setup.prior_info.delay.mean=58;
experiment_setup.prior_info.delay.std=15;
experiment_setup.prior_info.delay.n_grid=200;


load('l23_template_cell.mat');
temp = l23_average_shape; temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;
shape_template=l23_average_shape;


experiment_setup.prior_info.template_cell=struct;
experiment_setup.prior_info.template_cell.V_reset=-1e5;
experiment_setup.prior_info.template_cell.V_threshold=15;
experiment_setup.prior_info.template_cell.membrane_resistance = 0.02;
experiment_setup.prior_info.template_cell.cell_shape = shape_template;
experiment_setup.prior_info.template_cell.linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};

experiment_setup.prior_info.induced_intensity=struct;
experiment_setup.prior_info.induced_intensity.max_actual_stimulation=5;
experiment_setup.prior_info.induced_intensity.num_stim_grid=1000;
experiment_setup.prior_info.induced_intensity.linkfunc={@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};
experiment_setup.prior_info.induced_intensity.stim_scale = experiment_setup.prior_info.induced_intensity.num_stim_grid/experiment_setup.prior_info.induced_intensity.max_actual_stimulation;
experiment_setup.prior_info.induced_intensity.stim_grid = (1:experiment_setup.prior_info.induced_intensity.num_stim_grid)/...
experiment_setup.prior_info.induced_intensity.stim_scale;
% we should have a file with the values below precomputed... or an option
% to load from file if our priors haven't changed...
experiment_setup.prior_info.induced_intensity=precalculate_intensity(experiment_setup.prior_info.induced_intensity,...
experiment_setup.prior_info.template_cell,experiment_setup.prior_info.delay,experiment_setup.prior_info.current_template);
%         experiment_setup.prior_info.induced_intensity.intensity_grid
%         experiment_setup.prior_info.induced_intensity.probility_grid
%         experiment_setup.prior_info.induced_intensity.minimum_stim_threshold
%         experiment_setup.prior_info.induced_intensity.fire_stim_threshold


experiment_setup.neighbourhood_params=struct;
experiment_setup.neighbourhood_params.number=10;
experiment_setup.neighbourhood_params.height=15;
experiment_setup.neighbourhood_params.buffer_height=5; 



%
%
% %----------- Design parameters
% experiment_setup.design.num_groups = 3;
% experiment_setup.design.n_spots_per_trial = 3;
% experiment_setup.design.n_replicates=1;
% experiment_setup.design.K_undefined=8;
% experiment_setup.design.K_disconnected=12;
% experiment_setup.design.K_connected=4;
% experiment_setup.design.reps_undefined_single=8;
% experiment_setup.design.reps_disconnected_single=12;
% experiment_setup.design.reps_connected=4;
%
% experiment_setup.design.stim_loc_type = 1;
% experiment_setup.r1=5;experiment_setup.r2=10;
% experiment_setup.num_per_grid=12;
% experiment_setup.num_per_grid_dense = 16;
%
% experiment_setup.design.single_spot_threshold=9; % switch to single spot stimulation if there are fewer than N cells in this group
% experiment_setup.design.trial_max=20000;
% experiment_setup.design.disconnected_threshold = 0.2;
% experiment_setup.design.disconnected_confirm_threshold = 0.2;
%
%
% experiment_setup.design.connected_threshold = 0.5;
% experiment_setup.design.connected_confirm_threshold = 0.5;
%
% experiment_setup.design.n_MC_samples = 25;
%
%
% % Prior distribution
% experiment_setup.design.prior_pi0=0.8;
%
%
% % Initialize the variational family
% experiment_setup.design.gain_bound.up=0.03;
% experiment_setup.design.gain_bound.low=0.005;
%
% experiment_setup.design.var_pi_ini=0.01;% not used.
% experiment_setup.design.var_alpha_initial=1;
% experiment_setup.design.var_beta_initial=1.78;
% experiment_setup.design.var_alpha_gain_initial=...
%     log( (0.02 - experiment_setup.design.gain_bound.low)./(experiment_setup.design.gain_bound.up-0.02));
% experiment_setup.design.var_beta_gain_initial=0.5;
%
% % Initialize the parameters in the VI
%
% experiment_setup.design.C_threshold = 0.01;experiment_setup.design.maxit=1000;
% experiment_setup.design.S=200;experiment_setup.design.epsilon=0.01;experiment_setup.design.eta_logit=0;experiment_setup.design.eta_beta=0.05;
% experiment_setup.design.background_rt=experiment_setup.bg_rate*(experiment_setup.time.max_time - experiment_setup.time.min_time);
%
%
% experiment_setup.design.prob_weight=0;
%
% experiment_setup.design.lklh_func=@calculate_likelihood_sum_bernoulli;
% experiment_setup.design.stim_threshold = 10;


% experiment_setup.design.id_notconnected=false;
% experiment_setup.design.connected=true;
%
%  %  loc_to_cell_nuclei is from get_stim_locations
%
% experiment_setup.design.change_threshold=0.05;
% experiment_setup.design.do_connected_vi = 1;



% some experimental experiment_setup
% experiment_setup.exp.power_levels = '20 30 40 50 60 70'; % this should be a space delimited string
% experiment_setup.power_level=30:10:100;
% experiment_setup.exp.power_levels = mat2str(experiment_setup.power_level);
% experiment_setup.exp.power_levels = experiment_setup.exp.power_levels(2:end-1);% remove brackets
experiment_setup.exp.z_width = 20;
% experiment_setup.exp.z_depths = '10 30 50 70 90';% this should be a space delimited string
experiment_setup.exp.arbitrary_z = 0;

load('power-calibration.mat');
experiment_setup.exp.ratio_map = ratio_map;
experiment_setup.exp.pockels_lut = pockels_lut;
experiment_setup.exp.max_power_ref = max(experiment_setup.exp.pockels_lut(2,:));


experiment_setup.exp.max_spike_freq = .5; % don't revisit cells on average sooner than this in Hz
experiment_setup.exp.max_stim_freq = 20; % max frequency for system
experiment_setup.exp.foe_bounds = [-148 148; -148 148; 0 400];

% more experimental experiment_setup
experiment_setup.exp.sim_locs = 0;
experiment_setup.exp.stim_duration = .003; % length of laser on

experiment_setup.phase_base_file = 0;
if strcmp(experiment_setup.experiment_type,'experiment') || experiment_setup.sim.compute_phase_masks
    if experiment_setup.phase_base_file
        load('phase-mask-base.mat');
    else % load from base ws
        disk_grid_phase = evalin('base','disk_grid_phase');
        disk_grid_key = evalin('base','disk_grid_key');
        fine_spots_grid_phase = evalin('base','fine_spots_grid_phase');
        fine_spots_grid_key = evalin('base','fine_spots_grid_key');
        
        % 2p image alignment params
        image_zero_order_coord = evalin('base','image_zero_order_coord');
        image_um_per_px = evalin('base','image_um_per_px');
        stack_um_per_slice = evalin('base','stack_um_per_slice');
        
    end
    experiment_setup.disk_grid_phase = disk_grid_phase;
    experiment_setup.disk_grid_key = disk_grid_key;
    experiment_setup.fine_spots_grid_phase = fine_spots_grid_phase;
    experiment_setup.fine_spots_grid_key = fine_spots_grid_key;
    
    experiment_setup.image_zero_order_coord = image_zero_order_coord;
    experiment_setup.image_um_per_px = image_um_per_px;
    experiment_setup.stack_um_per_slice = stack_um_per_slice;
    
else
    experiment_setup.disk_grid_phase = [];
    experiment_setup.disk_grid_key = [];
    experiment_setup.fine_spots_grid_phase = [];
    experiment_setup.fine_spots_grid_key = [];
    
end

experiment_setup.exp.max_trials_per_sweep = 1250;
experiment_setup.exp.first_stim_time = 1.0; % in sec
experiment_setup.exp.filter_config = 'Femto Phasor';
experiment_setup.exp.sweep_time_padding = 5.0; % in sec
