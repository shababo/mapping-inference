%% Load the prior information
% run('./Simulation_preprocessing.m')

load('../Environments/l23_template_cell.mat');
temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;
shape_template=l23_average_shape;

load('../Environments/chrome-template-3ms.mat');
downsamp=1;time_max=300;
current_template=template(1:downsamp:time_max);


prior_info=struct;
    prior_info.template_cell=struct;
        prior_info.template_cell.V_reset=-1e5;
        prior_info.template_cell.V_threshold=15;
        prior_info.template_cell.membrane_resistance =0.02;
        prior_info.template_cell.cell_shape = shape_template;
    prior_info.PR_prior = struct;
        prior_info.PR_prior.type=1; % spike and logitnorm slab 
        prior_info.PR_prior.pi_logit=0;
        prior_info.PR_prior.alpha=0;
        prior_info.PR_prior.beta=1;
    prior_info.current_template = current_template;
    prior_info.gain_model=[];
    prior_info.delay=struct;
        prior_info.delay.delayed=true;    
        prior_info.delay.type='Gamma';
        prior_info.delay.mean=58;
        prior_info.delay.std=15;
        prior_info.delay.n_grid=200;


    %% Load the data set for cell locations
load('../Environments/6_3_s2c3_mpp_and_stim_data.mat')
cell_locations=cell_locs;
n_cell = size(cell_locations,1);
clear mpp stim_pow target_inds;

%% Simulation some more neurons:
switch spatial_density %
    case 1 % default, original cells
        n_extra=0;
    case 2 % simulate 50% more neurons
        n_extra=ceil(n_cell/2);
    case 3 % simulate twice as many neurons
        n_extra=ceil(n_cell);    
end
extra_locations = 2*pi*rand(n_extra,3); % only need the first two columns 
for i= 1:n_extra
    temp_cell=randsample(1:n_cell,1);
    extra_locations(i,:)=cell_locations(temp_cell,:)+ ...
        distance_neighbour*[sin(extra_locations(i,1))*sin(extra_locations(i,2)),cos(extra_locations(i,1))*sin(extra_locations(i,2)),cos(extra_locations(i,2))];
end

% for i=1:3
%     extra_locations(:,i)=extra_locations(:,i)*range(cell_locations(:,i))+min(cell_locations(:,i));
% end

cell_locations=[cell_locations;extra_locations];

n_cell = size(cell_locations,1);
%% Generate the cell parameters
background_rate=1e-4;

switch sparsity
    case 1 % very sparse
        connected_proportion=0.1;
    case 2 % normal
        connected_proportion=0.2;
    case 3 % dense
        connected_proportion=0.4;
end


funcs.invlink = @invlink_sig;%@(resp) log(1 + exp(resp));%@(x) exp(x);
switch cell_parameters_type
    case 1 % normal
        gamma_truth = (rand([n_cell 1])<connected_proportion).*(0.5+0.3*rand([n_cell 1]));
        gain_truth=0.02+(rand([n_cell 1])-0.5)*0.01;
        disp('normal gains and gammas')
    case 2 % normal gamma, extreme gains
        gain_bound_truth.up=0.03;gain_bound_truth.low=0.005;
        
        gamma_truth = (rand([n_cell 1])<connected_proportion).*(0.5+0.3*rand([n_cell 1]));
        n_connected=sum(gamma_truth>0);
        lognormal_temp=exp(normrnd(0,sqrt(1),[n_cell 1]));
        gain_truth= (gain_bound_truth.up-gain_bound_truth.low)*exp(-lognormal_temp)./(1+exp(-lognormal_temp))+...
            gain_bound_truth.low;
        gain_truth=max(gain_bound_truth.low+0.0005,gain_truth);
        disp('log-normal gains, normal gammas')
    case 3 % weak gamma, normal gains
        weak_gamma_proportion=0.5; % some gammas are small
        gamma_truth = (rand([n_cell 1])<connected_proportion).*(0.5+0.3*rand([n_cell 1]));
        n_connected=sum(gamma_truth>0);
        weak_index=randsample(find(gamma_truth>0), floor(weak_gamma_proportion*n_connected));
        gamma_truth(weak_index)=(0.15+0.1*rand([length(weak_index) 1]));
        gain_truth=0.02+(rand([n_cell 1])-0.5)*0.01;
        disp('normal gains, many weak gammas')
        
end


%%
cell_params=struct([]);
for i_cell = 1:n_cell
    cell_params(i_cell).v_thresh= 15;
    cell_params(i_cell).v_reset= -1e4;
    cell_params(i_cell).location=cell_locations(i_cell,:);
    cell_params(i_cell).g=0.02;
    cell_params(i_cell).optical_gain=gain_truth(i_cell);
    cell_params(i_cell).gamma=gamma_truth(i_cell);
    cell_params(i_cell).shape_truth=1; % for simulation
    cell_params(i_cell).shape=1;
    
end

%% Put the cells into ten non-overlapping groups by their z-coordinates

z_thresholds = quantile(cell_locations(:,3), (1:(n_planes))*0.1);
% alternatively, the thresholds can be determined by absolute depths
for i_cell = 1:size(cell_locations,1)
    cell_params(i_cell).group= sum(cell_locations(i_cell,3)>z_thresholds)+1;
end

%% Calculate the z-planes:
z_planes=zeros(n_planes,1);
for i_plane = 1:n_planes
    
    temp_loc=reshape([cell_params([cell_params(:).group]==i_plane ).location],3,[])';
    z_planes(i_plane)= mean(temp_loc(:,3)); 
end

related_cells_by_plane=cell([n_planes,1]);
temp_loc=reshape([cell_params(:).location],3,[])';
all_z=temp_loc(:,3);
for i_plane = 1:n_planes
    related_cells_by_plane{i_plane}= find( abs(all_z-z_planes(i_plane))<distance_plane  ); 
end
%% Select one plane as the primary plane
%
n_cell_this_plane=sum([cell_params(:).group]==this_plane);
target_cell_list=struct([]);
target_cell_list(1).primary=find([cell_params(:).group]==this_plane);

% include also the cells in neighbouring planes to account for their inputs
% neighbour_plane = [this_plane-1 this_plane+1];
% target_cell_list(1).secondary=find(ismember([cell_params(:).group],neighbour_plane));

% Include cells that are close to the chosen plane:

target_cell_list(1).secondary=setdiff(related_cells_by_plane{this_plane},target_cell_list(1).primary)';


%% For simulation
% related_cell_list=[target_cell_list.primary; target_cell_list.secondary];
% v_th_known_related= v_th_known(related_cell_list);
% v_reset_known_related=v_reset_known(related_cell_list);
% g_related = g_truth(related_cell_list);
% gamma_related = gamma_truth(related_cell_list);
% gain_related=gain_truth(related_cell_list);

                
%% Specify settings in this experiment 
% Some fields are filled during the experiment 

experiment_setup=struct; 
    experiment_setup.experiment_type='Simulation'; % Experiment; Simulation; Reproduction
   
    
%% Group specifications:

%% The undefined group:
    undefined_profile=struct;
        undefined_profile.group_ID=1;
        undefined_profile.group_type_ID='Undefined';
        undefined_profile.design_function=@design_undefined;
        undefined_profile.design_func_params=struct;
            undefined_profile.design_func_params.candidate_grid_params=struct;
                undefined_profile.design_func_params.candidate_grid_params.radius=[5 10];
                undefined_profile.design_func_params.candidate_grid_params.number=[12 16];
                undefined_profile.design_func_params.candidate_grid_params.grid_type='ring';
                % 2d ring or 3d sphere
            undefined_profile.design_func_params.trials_params=struct;
                undefined_profile.design_func_params.trials_params.replicates=1;
                undefined_profile.design_func_params.trials_params.spots_per_trial=4;
                undefined_profile.design_func_params.min_trials_per_cell=10;
                undefined_profile.design_func_params.min_trials_per_cell=10;
                undefined_profile.design_func_params.trials_params.power_level=30:10:100;
                undefined_profile.design_func_params.trials_params.stim_design='Optimal';
                undefined_profile.inference_params.MCsamples_for_posterior=50;
                % Random, Nuclei, or Optimal
                undefined_profile.design_func_params.trials_params.weighted_indicator=true;
                %   whether to conduct more trials on low PR cells
        undefined_profile.inference_function = @inference_undefined;
        undefined_profile.inference_params=struct;
            undefined_profile.inference_params.likelihood=@calculate_loglikelihood_bernoulli;
            undefined_profile.inference_params.maxit=1e4;
            undefined_profile.inference_params.MCsamples_for_gradient=50;
            undefined_profile.inference_params.convergence_threshold=1e-3;
            undefined_profile.inference_params.step_size=1;
            undefined_profile.inference_params.step_size_max=2;
            undefined_profile.inference_params.background_rate=1e-4;
            % this should be estimated in the experiment 
            undefined_profile.inference_params.MCsamples_for_posterior=50;
            undefined_profile.inference_params.recent_batches=2;
            undefined_profile.inference_params.bounds=struct;
                undefined_profile.inference_params.bounds.PR=[0.05 1];
                undefined_profile.inference_params.bounds.gain=[0.005 0.03];
         undefined_profile.regroup_function=@regroup_undefined;
         undefined_profile.regroup_func_params=struct;
            undefined_profile.regroup_func_params.connected_threshold=0.5;
            undefined_profile.regroup_func_params.disconnected_threshold=0.2;
            undefined_profile.regroup_func_params.quantile_prob=[0.05 0.95];
            undefined_profile.regroup_func_params.regroup_type='Quantiles';
            % Quantiles or NonzeroProb
            undefined_profile.regroup_func_params.singlespot_threshold=0.2;
            % certain proportion of cells 
            undefined_profile.regroup_func_params.change_threshold =0.05;
            

%% Neighbourhood initialization
neighbourhoods=struct;
neighbourhoods(1)=struct;
    neighbourhoods(1).neighbourhood_ID=1;
    neighbourhoods(1).neurons=struct;
        neighbourhoods(1).neurons(i)
            neighbourhoods(1).neurons(i).cell_ID
            neighbourhoods(1).neurons(i).group_type_ID
            neighbourhoods(1).neurons(i).location
            neighbourhoods(1).neurons(i).V_reset
            neighbourhoods(1).neurons(i).V_thresh
            neighbourhoods(1).neurons(i).membrane_resistance
            neighbourhoods(1).neurons(i).PR_params=struct;
                neighbourhoods(1).neurons(i).PR_params(batch_ID).pi=0;
                neighbourhoods(1).neurons(i).PR_params(batch_ID).alpha=0;
                neighbourhoods(1).neurons(i).PR_params(batch_ID).beta=0;
                neighbourhoods(1).neurons(i).PR_params(batch_ID).type='SpikedLogitNormal';
                neighbourhoods(1).neurons(i).PR_params(batch_ID).mean=0;
                neighbourhoods(1).neurons(i).PR_params(batch_ID).variance=0;
                neighbourhoods(1).neurons(i).PR_params(batch_ID).upper_quantile=0;
                neighbourhoods(1).neurons(i).PR_params(batch_ID).lower_quantile=0;
                neighbourhoods(1).neurons(i).PR_params(batch_ID).nonzero_prob=0;
        neighbourhoods(1).neurons(i).gain_params;
        neighbourhoods(1).neurons(i).primary_indicator=true;
        % whether it is a cell in this neighbourhood, or a nearby cell
        neighbourhoods(1).neurons(i).stim_locations=struct;
            neighbourhoods(1).neurons(i).stim_locations.undefined=struct;
                neighbourhoods(1).neurons(i).stim_locations.undefined.grid=[];
                neighbourhoods(1).neurons(i).stim_locations.undefined.effect=[];
        neighbourhoods(1).neurons(i).truth=struct; % parameters for simulation
        neighbourhoods(1).computing_time=struct;
            neighbourhoods(1).computing_time(batch_ID).undefined=0;
            neighbourhoods(1).computing_time(batch_ID).connected=0;
%% Designing experiment
%-------------------------------------------------%
%% Parameters in the design stage
n_replicates=1; % number of replicates for each trial
n_spots_per_trial = 4;trial_max=2000;
K_connected=10;K_undefined=20;

switch assignment_type
    case 1
        disconnected_threshold = 0.1;connected_threshold = 0.4;
        quantiles_prob=[0.1;0.9];
    case 2
        disconnected_threshold = 0.05;connected_threshold = 0.95;
        gamma_bound.up=1;gamma_bound.low=0.05;
        quantiles_prob=[0.1;0.90];
    otherwise
end
        
% if the postsynaptic cell has good connection, switch to 
disconnected_threshold_good = 0.2;connected_threshold_good = 0.5;
quantiles_prob_good=[0.1;0.90];

multiplier_connected=1;
single_spot_threshold=12; %12*10/4=30 (only 30 trials on these cells)
% switch to single spot stimulation (this can be a function of n_spots_per_trial
multispot_change_threshold=0.05; % switch if there aren't enough changes 

% Define tuning parameters in the VI
maxit=1000;
S=50;epsilon=0.01;eta=1;eta_max=2;
background_rt=background_rate*time_max; % raw probability of firing within a trial
n_MC_samples=50;


num_trials_per_batch=200;
% tuning parameter: prior_info_type;
run('./Simulation_tuning.m')
spike_indicator=false;
%%

disconnected_params=struct([]);
disconnected_params(1).prop=0.8;
disconnected_params(1).alpha=-4;
disconnected_params(1).beta=0;

regular_params=struct([]);
regular_params(1).prop=0.2;
regular_params(1).alpha=0.6;
regular_params(1).beta=0;

som_params=struct([]);
som_params(1).prop=0.2;
som_params(1).alpha=-1.4;
som_params(1).beta=0;




% tuning parameters: model_type, design_type, assignment_type
run('./Simulation_experiment.m')

