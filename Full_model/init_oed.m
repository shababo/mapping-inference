function params = init_oed

% parameters


%----------- Delay parameters
params.delay.type=2; %1: normal; 2: gamma
params.delay.mean=58; delay_params.std=15;
params.delay.delayed=true; params.delay.n_grid=200;

params.bg_rate = 1e-4;

%----------- Load the current template
load('./Environments/chrome-template-3ms.mat');
params.time.downsamp=1;params.time.max_time=300;
params.current_template=template(1:downsamp:max_time);
params.t_vect= 1:1:params.time.max_time;

%----------- Build template cell
params.template_cell.gain_template=0.02;
params.template_cell.g=0.02;params.template_cell.v_th_known=15;
params.template_cell.v_reset_known = 1e-4;
params.template_cell.stim_unique = (1:1000)/10;
params.template_cell.linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};
load('./Environments/l23_template_cell.mat');
temp=l23_average_shape;temp_max = max(max(max(temp)));
params.template_cell.shape_template=temp/temp_max;
params.stim_unique = (1:1000)/10;

[params.template_cell.prob_trace]=get_firing_probability(...
    params.template_cell.linkfunc,params.current_template,params.stim_unique,params.template_cell,params.delay);

% Calculate the firing intensity
params.stim_scale=200;
params.stim_grid = (1:1000)/stim_scale;
% cell_params.g=0.02;
[params.template_cell.prob_trace_full,params.template_cell.v_trace_full] = get_first_spike_intensity(...
    params.template_cell.linkfunc,...
    params.current_template,params.stim_grid,params.template_cell,params.delay);

%
params.eff_stim_threshold=stim_grid(min(find(sum(params.template_cell.prob_trace_full,2)>1e-1)));

%----------- Design parameters
params.design.num_groups = 3;
params.design.n_spots_per_trial = 3;
params.design.n_replicates=1; % conduct two replicates for each trial
params.design.K_undefined=6; % each cell appears approximately 10*2 times
params.design.K_disconnected=6; % each cell appears approximately 10*2 times
params.design.K_connected=20; % each cell appears approximately 10*2 times

params.design.single_spot_threshold=15; % switch to single spot stimulation if there are fewer than N cells in this group
params.design.trial_max=2000;
params.design.disconnected_threshold = 0.2;
params.design.disconnected_confirm_threshold = 0.2;


params.design.connected_threshold = 0.5;
params.design.connected_confirm_threshold = 0.5;




% Prior distribution
params.design.prior_pi0=0.8;


% Initialize the variational family
params.design.var_pi_ini=0.01;% not used.
params.design.var_alpha_initial=1;params.design.var_beta_initial=1.78;
params.design.var_alpha_gain_initial=1;params.design.var_beta_gain_initial=1.78;

% Initialize the parameters in the VI
params.design.C_threshold = 0.01;params.design.maxit=1000;
params.design.S=200;params.design.epsilon=0.01;params.design.eta_logit=0;params.design.eta_beta=0.05;
params.design.background_rt=params.bg_rate*params.time.time_max;

params.design.gamma_estimates = 0.5*ones(n_cell_this_plane,1);% for drawing samples...


params.design.prob_weight=0;

params.design.lklh_func=@calculate_likelihood_sum_bernoulli;
params.design.stim_threshold = 10;
params.design.gain_bound.up=0.03;
params.design.gain_bound.low=0.01;

params.design.id_notconnected=false;
params.design.connected=true;

 %  loc_to_cell_nuclei is from get_stim_locations 

params.design.change_threshold=0.05;

% some experimental params
params.exp.power_levels = '50'; % this should be a space delimited string
params.exp.z_width = 30;
params.exp.ratio_map = evalin('base','ratio_map');
params.exp.pockels_lut = evalin('base','pockels_lut');
params.exp.max_ratio_ref = max(params.exp.pockels_lut)/...
    ((params.exp.power_level+15)*params.design.n_spots_per_trial);
params.exp.max_spike_freq = .5; % don't revisit cells on average sooner than this in Hz

