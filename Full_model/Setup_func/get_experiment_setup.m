
function experiment_setup = get_experiment_setup(varargin)

% read this from argument
group_names=cell(3 1);
group_names{1}='undefined';
group_names{2}='connected';
group_names{3}='disconnected';



%experiment_setup.is_sim

% parameters
if ~isempty(varargin) && ~isempty(varargin{1})
    load_map = varargin{1};
else
    load_map = 0;
end


experiment_setup.savedir = 'C:\data\Shababo';
experiment_setup.root = '/media/shababo/data';

clock_array = clock;
experiment_setup.map_id = [num2str(clock_array(2)) '_' num2str(clock_array(3)) ...
    '_' num2str(clock_array(4)) ...
    '_' num2str(clock_array(5))];
experiment_setup.exp_id = experiment_setup.map_id;
experiment_setup.fullsavefile = fullfile(experiment_setup.savedir,[experiment_setup.map_id '_data.mat']);

%----------- Delay parameters
experiment_setup.delay.type=2; %1: normal; 2: gamma
experiment_setup.delay.mean=58; experiment_setup.delay.std=15;
experiment_setup.delay.delayed=true; experiment_setup.delay.n_grid=200;

experiment_setup.bg_rate = 1e-4;

%----------- Load the current template

experiment_setup.experiment_type='Simulation'; % Experiment; Simulation; Reproduction
experiment_setup.prior_info=struct;

experiment_setup.prior_info.PR_prior = struct;
experiment_setup.prior_info.PR_prior.type=1; % spike and logitnorm slab
experiment_setup.prior_info.PR_prior.pi_logit=0;
experiment_setup.prior_info.PR_prior.alpha=0;
experiment_setup.prior_info.PR_prior.beta=1;

load('../Environments/l23_template_cell.mat');
temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;
shape_template=l23_average_shape;

load('../Environments/chrome-template-3ms.mat');
downsamp=1;time_max=300;
current_template=template(1:downsamp:time_max);


experiment_setup.prior_info.current_template = current_template;
experiment_setup.prior_info.gain_model=[];
experiment_setup.prior_info.delay=struct;
experiment_setup.prior_info.delay.delayed=true;
experiment_setup.prior_info.delay.type='gamma';
experiment_setup.prior_info.delay.mean=58;
experiment_setup.prior_info.delay.std=15;
experiment_setup.prior_info.delay.n_grid=200;


load('l23_template_cell.mat');
temp=l23_average_shape;temp_max = max(max(max(temp)));

experiment_setup.prior_info.template_cell=struct;
experiment_setup.prior_info.template_cell.V_reset=-1e5;
experiment_setup.prior_info.template_cell.V_threshold=15;
experiment_setup.prior_info.template_cell.membrane_resistance =0.02;
experiment_setup.prior_info.template_cell.cell_shape = shape_template;
experiment_setup.prior_info.template_cell.linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};

experiment_setup.prior_info.induced_intensity=struct;
experiment_setup.prior_info.induced_intensity.max_actual_stimulation=5;
experiment_setup.prior_info.induced_intensity.num_stim_grid=1000;
experiment_setup.prior_info.induced_intensity.linkfunc={@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};
experiment_setup.prior_info.induced_intensity.stim_scale = experiment_setup.prior_info.induced_intensity.num_stim_grid/experiment_setup.prior_info.induced_intensity.max_actual_stimulation;
experiment_setup.prior_info.induced_intensity.stim_grid = (1:num_stim_grid)/stim_scale;
experiment_setup.prior_info.induced_intensity=precalculate_intensity(experiment_setup.prior_info.induced_intensity,...
    experiment_setup.prior_info.template_cell,experiment_setup.prior_info.delay,experiment_setup.prior_info.current_template);
%         experiment_setup.prior_info.induced_intensity.intensity_grid
%         experiment_setup.prior_info.induced_intensity.probility_grid
%         experiment_setup.prior_info.induced_intensity.minimum_stim_threshold
%         experiment_setup.prior_info.induced_intensity.fire_stim_threshold


experiment_setup.neighbourhood_params=struct;
experiment_setup.neighbourhood_params.number=10;
experiment_setup.neighbourhood_params.height=20;
experiment_setup.neighbourhood_params.buffer_height=5;
experiment_setup.neighbourhood_params.z_thresholds;


% need a way to call these functions based on group names
experiment_setup.groups.undefined=get_undefined();
experiment_setup.groups.connected=get_connected();
experiment_setup.groups.disconnected=get_disconnected();
experiment_setup.groups.alive=get_alive();
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

% for std thresh experiments
experiment_setup.design.std_thresh = [0 .50];
experiment_setup.design.min_targs = 10;

% some experimental experiment_setup
% experiment_setup.exp.power_levels = '20 30 40 50 60 70'; % this should be a space delimited string
experiment_setup.exp.power_levels = mat2str(experiment_setup.power_level);
experiment_setup.exp.power_levels = experiment_setup.exp.power_levels(2:end-1);
experiment_setup.exp.z_width = 20;
experiment_setup.exp.z_depths = '10 30 50 70 90';% this should be a space delimited string
experiment_setup.exp.arbitrary_z = 0;

if load_map
    experiment_setup.exp.ratio_map = evalin('base','ratio_map');
    experiment_setup.exp.pockels_lut = evalin('base','pockels_lut');
    
    experiment_setup.exp.max_power_ref = max(experiment_setup.exp.pockels_lut(2,:));
end

experiment_setup.exp.max_spike_freq = .5; % don't revisit cells on average sooner than this in Hz
experiment_setup.exp.max_stim_freq = 15;
experiment_setup.exp.foe_bounds = [-148 148; -148 149];

% more experimental experiment_setup
experiment_setup.exp.sim_locs = 0;

experiment_setup.exp.rand_order = 1;
% set(handles.num_repeats,'String',num2str(10));
experiment_setup.exp.duration = .003; % length of laser on


% if filename given
% load filename
% for every field in params from file
% replace default with file version


