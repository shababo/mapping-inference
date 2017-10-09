%% Preprocessing
%--------------------------------------------------%
%% Estimate the first spike probability given a template cell:
%----------- Delay parameters
delay_params.type=2; %1: normal; 2: gamma
delay_params.mean=58; delay_params.std=15;
delay_params.delayed=true; delay_params.n_grid=200;

cell_params.gain_template = 1; % for calculation
cell_params.g=0.02;cell_params.v_th_known=15;
linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};

gain_template=0.02;
stim_scale=4/gain_template;
stim_grid = (1:1000)/stim_scale;
stim_unique=(1:1000)/stim_scale/gain_template;
[prob_trace_full,v_trace_full] = get_first_spike_intensity(...
    linkfunc,...
    current_template,stim_grid,cell_params,delay_params);
prob_trace=sum(prob_trace_full,2);
eff_stim_threshold=stim_grid(min(find(prob_trace>1e-1)));
spiked_stim_threshold = stim_grid(min(find(prob_trace>8e-1)));
%% Select the stimulation locations
% Load the shape template
load('../Environments/l23_template_cell.mat');
load('../Environments/l23_cells_for_sim.mat');
%l23_average_shape
temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;shape_template=l23_average_shape;
r1=5;r2=10;r3=15;num_per_grid=12;
num_per_grid_dense=16;
stim_threshold=eff_stim_threshold/gain_template;
grid_type = 1; % 1: circles; 2: lines;


[~,~,~,~,...
    target_locations_all,cell_neighbours,...
    target_locations_nuclei, power_nuclei,pi_target_nuclei, loc_to_cell_nuclei] = ...
    get_stim_locations(...
    target_cell_list,cell_locations,power_level,...
    r1,r2,r3,num_per_grid,num_per_grid_dense,shape_template,...
    stim_unique,prob_trace,stim_threshold,...
    grid_type);

%% Overwrite the selected locations 
target_locations_selected=cell_locations(target_cell_list.primary,:);
target_locations_selected(:,3)=mean(target_locations_selected(:,3));


cell_params.locations =  cell_locations(related_cell_list,:);
cell_params.shape_gain = ones(length(related_cell_list),1);
cell_template = struct();
cell_template.shape= shape_template;
[pi_target_selected, inner_normalized_products] = get_weights_v2(cell_params, ...
    cell_template,target_locations_selected);
loc_to_cell = 1:size(target_locations_selected,1);

%% Calculate the initial power levels for cells using prior information of the gain variables 
stimulation_required =  spiked_stim_threshold./((1./(1+exp(-variational_params_path.alpha_gain(1:n_cell_this_plane) )))*(gain_bound.up-gain_bound.low)+gain_bound.low);

% power_selected=100*ones(length(power_selected),1);
% Initialize as zeros:
%power_selected=stimulation_required./diag(pi_target_selected);
power_selected=0*ones(length(power_selected),1);

power_sd=10*ones(length(power_selected),1); % this should be related to the posterior distribution
%% End of preprocessing
%-------------------------------------------------%
