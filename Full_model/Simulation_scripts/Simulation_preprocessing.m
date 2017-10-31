%% Preprocessing
%--------------------------------------------------%
%%  Load the current template
load('../Environments/chrome-template-3ms.mat');
downsamp=1;time_max=300;power_level = 30:10:100;
current_template=template(1:downsamp:time_max);
t_vect= 1:1:time_max;
%% Specify delay distributions
delay_params.type=2; %1: normal; 2: gamma
delay_params.mean=58; delay_params.std=15;
delay_params.delayed=true; delay_params.n_grid=200;

% bounds of the gamma:
gain_bound.up=0.03;gain_bound.low=0.003;
gamma_bound.up=1;gamma_bound.low=0.01;

%% Pre-calculate the first spike intensity as a function of actual stimulation 
template_cell.optical_gain = 1; % for calculation
template_cell.g=0.02;template_cell.v_thresh=15;

max_actual_stimulation=5;
num_stim_grid=1000;
linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};

gain_template=0.02;
stim_scale=num_stim_grid/max_actual_stimulation;
stim_grid = (1:num_stim_grid)/stim_scale; % grid of actual stimulation 

% stim_unique=(1:1000)/stim_scale/gain_template;
[prob_trace_full,~] = get_first_spike_intensity(...
    linkfunc,current_template,stim_grid,template_cell,delay_params);
prob_trace=sum(prob_trace_full,2);

% Minimum actual stimulation in order for a cell to fire 
minimum_stim_threshold=stim_grid(min(find(prob_trace>0.01))); 

% The threshold of actual stimulation when a cell fires with high
% probability 
fire_stim_threshold=stim_grid(min(find(prob_trace>0.99)));
stim_threshold = minimum_stim_threshold/gain_bound.up;

%% Pre-calculate the candidate stimulation locations 
load('../Environments/l23_template_cell.mat');
% load('../Environments/l23_cells_for_sim.mat');

temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;
shape_template=l23_average_shape;

grid_params=struct([]);
grid_params(1).radius=[5];
grid_params(1).number=[12];
grid_params(2).radius=[10];
grid_params(2).number=[16];


[pi_target, inner_normalized_products,target_locations,loc_to_cell,...
    target_locations_nuclei,pi_target_nuclei, loc_to_cell_nuclei]= ...
   get_stim_locations(target_cell_list,cell_params,grid_params,shape_template);


%% For simulation:
sim_indicator=true;
[pi_target_truth, ~,~,~,...
    ~,pi_target_nuclei_truth,~]= ...
   get_stim_locations(target_cell_list,cell_params,grid_params,shape_template,[],[],sim_indicator);

%% End of preprocessing
%-------------------------------------------------%
