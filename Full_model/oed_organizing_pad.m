%% data/input
cell_locations
power_level = [50 75 100];

%% parameters

params.n_cell = size(cell_locations,1);

%----------- Delay parameters
params.delay.type=2; %1: normal; 2: gamma
params.delay.mean=58; delay_params.std=15;
params.delay.delayed=true; params.delay.n_grid=200;

%----------- Load the current template
load('./Environments/chrome-template-3ms.mat');
downsamp=1;max_time=300;power_level = [50 75 100];
num_power_level=length(power_level);
current_template=template(1:downsamp:max_time);
t_vect= 1:1:max_time;

%----------- Build template cell
gain_template=0.02;
template_cell.g=0.02;template_cell.v_th_known=15;template_cell.gain_template = gain_template;
stim_unique = (1:1000)/10;
linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};
load('./Environments/l23_template_cell.mat');
temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max; template_cell.shape_template=l23_average_shape;

[prob_trace]=get_firing_probability(...
    linkfunc,current_template,stim_unique,template_cell,delay_params);


%%