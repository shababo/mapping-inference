function data = init_oed

% parameters

params.n_cell = size(cell_locations,1);

%----------- Delay parameters
params.delay.type=2; %1: normal; 2: gamma
params.delay.mean=58; delay_params.std=15;
params.delay.delayed=true; params.delay.n_grid=200;

%----------- Load the current template
load('./Environments/chrome-template-3ms.mat');
params.time.downsamp=1;params.timemax_time=300;
params.current_template=template(1:downsamp:max_time);
params.t_vect= 1:1:max_time;

%----------- Build template cell
params.template_cell.gain_template=0.02;
params.template_cell.g=0.02;params.template_cell.v_th_known=15;
params.template_cell.stim_unique = (1:1000)/10;
params.template_cell.linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};
load('./Environments/l23_template_cell.mat');
temp=l23_average_shape;temp_max = max(max(max(temp)));
params.template_cell.shape_template=temp/temp_max;

[params.template_cell.prob_trace]=get_firing_probability(...
    linkfunc,current_template,stim_unique,params.template_cell,delay_params);

% Calculate the firing intensity
params.stim_scale=200;
params.stim_grid = (1:1000)/stim_scale;
% cell_params.g=0.02;
[params.template_cell.prob_trace_full,params.template_cell.v_trace_full] = get_first_spike_intensity(...
    linkfunc,...
    current_template,stim_grid,params.template_cell,delay_params);

%
params.eff_stim_threshold=stim_grid(min(find(sum(prob_trace_full,2)>1e-1)));

data.params = params;