
function params = init_oed(varargin)

% parameters
if ~isempty(varargin) && ~isempty(varargin{1})
    load_map = varargin{1};
else
    load_map = 0;
end


params.savedir = 'C:\data\Shababo';

clock_array = clock;
params.map_id = [num2str(clock_array(2)) '_' num2str(clock_array(3)) ...
    '_' num2str(clock_array(4)) ...
    '_' num2str(clock_array(5))];
params.fullsavefile = fullfile(params.savedir,[params.map_id '_data.mat']);

%----------- Delay parameters
params.delay.type=2; %1: normal; 2: gamma
params.delay.mean=58; params.delay.std=15;
params.delay.delayed=true; params.delay.n_grid=200;

params.bg_rate = 1e-4;

%----------- Load the current template
load('chrome-template-3ms.mat');
params.time.downsamp=1;params.time.max_time=300;params.time.min_time = 45;
params.current_template=template(1:params.time.downsamp:params.time.max_time);
params.t_vect= 1:1:params.time.max_time;

params.power_level = [25 50 75 100];
params.num_power_level=length(params.power_level);

%----------- Build template cell
params.template_cell.gain_template=0.02;
params.template_cell.g=0.02;params.template_cell.v_th_known=15;
params.template_cell.v_reset_known = 1e-4;
params.template_cell.stim_unique = (1:1000)/10;
params.template_cell.linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};
load('l23_template_cell.mat');
temp=l23_average_shape;temp_max = max(max(max(temp)));
params.template_cell.shape_template=temp/temp_max;
params.stim_scale=4/params.template_cell.gain_template;
params.stim_grid = (1:1000)/params.stim_scale;
params.stim_unique = (1:1000)/params.stim_scale/params.template_cell.gain_template;

% [params.template_cell.prob_trace]=get_firing_probability(...
%     params.template_cell.linkfunc,params.current_template,params.stim_unique,params.template_cell,params.delay);

% Calculate the firing intensity


% cell_params.g=0.02;
[params.template_cell.prob_trace_full,params.template_cell.v_trace_full] = get_first_spike_intensity(...
    params.template_cell.linkfunc,...
    params.current_template,params.stim_grid,params.template_cell,params.delay);
params.template_cell.prob_trace=sum(params.template_cell.prob_trace_full,2);
%
% params.eff_stim_threshold=params.stim_grid(min(find(sum(params.template_cell.prob_trace_full,2)>1e-1)));
params.eff_stim_threshold=params.stim_grid(min(find(params.template_cell.prob_trace>0.01)));
params.fire_stim_threshold=params.stim_grid(min(find(params.template_cell.prob_trace>0.99)));

%----------- Design parameters
params.design.num_groups = 3;
params.design.n_spots_per_trial = 3;
params.design.n_replicates=1; % conduct two replicates for each trial
params.design.K_undefined=8; % each cell appears approximately 10*2 times
params.design.K_disconnected=8; % each cell appears approximately 10*2 times
params.design.K_connected=4; % each cell appears approximately 10*2 times
params.design.reps_undefined_single=8;
params.design.reps_disconnected_single=8;
params.design.reps_connected=4;

params.design.stim_loc_type = 1;

params.design.single_spot_threshold=9; % switch to single spot stimulation if there are fewer than N cells in this group
params.design.trial_max=20000;
params.design.disconnected_threshold = 0.2;
params.design.disconnected_confirm_threshold = 0.2;


params.design.connected_threshold = 0.5;
params.design.connected_confirm_threshold = 0.5;

params.design.n_MC_samples = 25;


% Prior distribution
params.design.prior_pi0=0.8;


% Initialize the variational family
params.design.var_pi_ini=0.01;% not used.
params.design.var_alpha_initial=1;params.design.var_beta_initial=1.78;
params.design.var_alpha_gain_initial=1;params.design.var_beta_gain_initial=1.78;

% Initialize the parameters in the VI
params.design.C_threshold = 0.01;params.design.maxit=1000;
params.design.S=200;params.design.epsilon=0.01;params.design.eta_logit=0;params.design.eta_beta=0.05;
params.design.background_rt=params.bg_rate*(params.time.max_time - params.time.min_time);


params.design.prob_weight=0;

params.design.lklh_func=@calculate_likelihood_sum_bernoulli;
params.design.stim_threshold = 10;
params.design.gain_bound.up=0.03;
params.design.gain_bound.low=0.01;

params.design.id_notconnected=false;
params.design.connected=true;

 %  loc_to_cell_nuclei is from get_stim_locations 

params.design.change_threshold=0.01;
params.design.do_connected_vi = 1;

% for std thresh experiments
params.design.std_thresh = [0 .50];
params.design.min_targs = 10;

% some experimental params
params.exp.power_levels = '50'; % this should be a space delimited string
params.exp.z_width = 20;
params.exp.z_depths = '30 50';% this should be a space delimited string
params.exp.arbitrary_z = 0;

if load_map
    params.exp.ratio_map = evalin('base','ratio_map');
    params.exp.pockels_lut = evalin('base','pockels_lut');
    
    params.exp.max_power_ref = max(params.exp.pockels_lut(2,:));
end

params.exp.max_spike_freq = .5; % don't revisit cells on average sooner than this in Hz
params.exp.max_stim_freq = 15;
params.exp.foe_bounds = [-148 148; -148 148];


