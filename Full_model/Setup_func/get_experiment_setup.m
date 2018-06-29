function experiment_setup = get_experiment_setup(varargin)

% to-dos:
% 1) load in the prior distributions from existing files via initialize
% parameters
% 2) plug in the gain ~ fluor + depth model
% experiment_setup.prior_info.gain_model

%experiment_setup.is_sim

% parameters

if ~isempty(varargin) && ~isempty(varargin{1})
    param_preset_ID = varargin{1};
else
    param_preset_ID = 'millennium-falcon-exp-sim';
end

% experiment; simulation; reproduction

 % experiment; simulation; reproduction

switch param_preset_ID
    case 'szchen-sim'
        
        experiment_setup.prior_root = 'C:/Users/Shizhe/Documents/Mapping_data/Data/';
        experiment_setup.exp_root = 'C:/Users/Shizhe/Documents/Mapping_data/Data/';
        experiment_setup.analysis_root = 'C:/Users/Shizhe/Documents/Mapping_data/tmp/';
        experiment_setup.experiment_type='simulation';
        experiment_setup.single_patch_path='C:/Users/Shizhe/Documents/Mapping_data/Data/more_cells.mat';
    case 'szchen-rep'
        experiment_setup.exp_root = 'C:/Users/Shizhe/Documents/Mapping_data/Data/';
        experiment_setup.analysis_root = 'C:/Users/Shizhe/Documents/Mapping_data/Data/NewData/';
        experiment_setup.experiment_type='reproduction';
        experiment_setup.single_patch_path='C:/Users/Shizhe/Documents/Mapping_data/Data/more_cells.mat';
        experiment_setup.rep.file_name=[experiment_setup.exp_root '/12_14_16_49/12_14_16_49_data.mat'];
    case 'szchen-sim-hab'
        experiment_setup.exp_root = '../Data/';
        experiment_setup.analysis_root = '../Data/AnalysisRoot/';
        experiment_setup.experiment_type='simulation';
        experiment_setup.single_patch_path='../Data/Environments/more_cells.mat';
    case 'szchen-rep-hab1'
        experiment_setup.exp_root = '../Data/';
        experiment_setup.analysis_root = '../Data/tmp';
        experiment_setup.experiment_type='reproduction';
        experiment_setup.single_patch_path='../Data/Environments/more_cells.mat';
        experiment_setup.rep.file_name=[experiment_setup.exp_root '12_14_16_49/12_14_16_49_data.mat'];
    case 'szchen-rep-hab2'
        experiment_setup.exp_root = '../Data/';
        experiment_setup.analysis_root = '../Data/tmp';
        experiment_setup.experiment_type='reproduction';
        experiment_setup.single_patch_path='../Data/Environments/more_cells.mat';
        experiment_setup.rep.file_name=[experiment_setup.exp_root '12_13_17_0/12_13_17_0_data.mat'];
    case 'millennium-falcon-mapping'
        experiment_setup.experiment_type='experiment';
        experiment_setup.exp_root = 'C:\data\Shababo\';
        experiment_setup.analysis_root = '/media/shababo/data/'; % make sure to add ending slash
        experiment_setup.ephys_mapped_drive = '/home/shababo/slidebook-comp/';
        experiment_setup.phase_mask_dir = 'W:\';
        experiment_setup.single_patch_path=[];
    case 'millennium-falcon-exp-sim'
        experiment_setup.experiment_type='experiment';
        experiment_setup.exp_root = 'C:\data\Shababo\';
        experiment_setup.analysis_root = '/media/shababo/data/'; % make sure to add ending slash
        experiment_setup.ephys_mapped_drive = '/home/shababo/slidebook-comp/';
        experiment_setup.phase_mask_dir = 'W:\';
        experiment_setup.single_patch_path=[];
    case 'turing-rep'
        experiment_setup.experiment_type='reproduction';
        experiment_setup.exp_root = '/media/shababo/data/';
        experiment_setup.analysis_root = '/media/shababo/data/'; % make sure to add ending slash
        experiment_setup.ephys_mapped_drive = '/home/shababo/slidebook-comp/';
        experiment_setup.phase_mask_dir = '/home/shababo/slidebook-comp/';
        experiment_setup.single_patch_path=[];
    case 'shababo-mbp-sim'
        experiment_setup.exp_root = '/Users/shababo/projects/mapping/data/sim_tmp/';
        experiment_setup.analysis_root = '/Users/shababo/projects/mapping/data/sim_tmp/';
        experiment_setup.single_patch_path=[];
end

clock_array = clock;
experiment_setup.exp_id = [num2str(clock_array(2)) '_' num2str(clock_array(3)) ...
    '_' num2str(clock_array(4)) ...
    '_' num2str(clock_array(5))];
experiment_setup.exp.fullsavefile = fullfile(experiment_setup.exp_root,[experiment_setup.exp_id '_data.mat']);

% read this from argument
group_names=cell([5 1]);
group_names{1}='undefined';
group_names{2}='connected';
group_names{3}='alive';
group_names{4}='disconnected';
group_names{5}='secondary';
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
switch  experiment_setup.experiment_type
    case {'simulation','experiment'}
        experiment_setup.sim = get_simulation_setup;
    case 'reproduction'
        experiment_setup.rep = get_reproduction_setup;
end

switch param_preset_ID
    case {'szchen-sim','szchen-sim-hab'}
        experiment_setup.sim.compute_phase_masks=0;
        experiment_setup.sim.use_power_calib =0;
    case {'szchen-rep','turing-rep'}
        experiment_setup.rep.compute_phase_masks=0;
        experiment_setup.rep.use_power_calib =0;
        experiment_setup.rep.file_name=[experiment_setup.analysis_root '12_14_16_49/12_14_16_49_data.mat'];
    case 'szchen-rep-hab1'
        experiment_setup.rep.file_name=[experiment_setup.exp_root '12_14_16_49/12_14_16_49_data.mat'];
    case 'szchen-rep-hab2'
        experiment_setup.rep.file_name=[experiment_setup.exp_root '12_13_17_0/12_13_17_0_data.mat'];
        
    case 'millennium-falcon-mapping'
        experiment_setup.sim.compute_phase_masks=1;
end



if strcmp(experiment_setup.experiment_type,'reproduction')
    if experiment_setup.rep.use_same_setting
        temp_setup=experiment_setup; % to preserve paths, etc.
        % NOTE: this will overwrite the experiment_setup
        load(temp_setup.rep.file_name)
        % Rewrite the path info on experiment_setup
        experiment_setup.experiment_type=temp_setup.experiment_type;
        experiment_setup.exp_root=temp_setup.exp_root;
        experiment_setup.analysis_root=temp_setup.analysis_root;
        experiment_setup.groups=temp_setup.groups;
       experiment_setup.rep=temp_setup.rep;
        if isfield(temp_setup,'sim')
            experiment_setup.sim=temp_setup.sim;
        end
        experiment_setup.records=struct;
        experiment_setup.records.queries=experiment_queries;
        experiment_setup.records.neighbourhoods=neighbourhoods;
        experiment_setup.reproduced=struct;
        % INITIALIZE WITH RECORDS TO AVOID DISSIMILAR STRUCTURES
        experiment_setup.reproduced.queries=experiment_queries;
        experiment_setup.reproduced.neighbourhoods=neighbourhoods;
    end
else
    
    experiment_setup.prior_info=struct;
    % Set up parameters for the prior distributions:
    % Note: neuron-specific prior & initial values are specified in
    % run_mapping_experiment
    experiment_setup.prior_info.prior_parameters=initialize_parameters();
    
    % load the GP parameters: 
    experiment_setup.prior_info.GP_params=initialize_GP(experiment_setup.prior_root);
    
  

    %----------- Load the current template
    load([ experiment_setup.prior_root 'chrome-template-3ms.mat']);
    
    
    experiment_setup.trials.downsamp = 1;
    experiment_setup.trials.Fs = 20000;
    experiment_setup.trials.min_time = 50;
    experiment_setup.trials.max_time = 500;
    experiment_setup.trials.max_time_sec = experiment_setup.trials.max_time/experiment_setup.trials.Fs;
    current_template=template(1:experiment_setup.trials.downsamp:experiment_setup.trials.max_time);
    
    experiment_setup.prior_info.current_template = current_template;
    experiment_setup.prior_info.gain_model=[];
%     experiment_setup.prior_info.delay=struct;
%     experiment_setup.prior_info.delay.delayed=true;
%     experiment_setup.prior_info.delay.type='gamma';
%     experiment_setup.prior_info.delay.mean=58;
%     experiment_setup.prior_info.delay.std=20;
%     experiment_setup.prior_info.delay.n_grid=200;
    
    load([ experiment_setup.prior_root 'l23_template_cell.mat']);
    temp = l23_average_shape; temp_max = max(max(max(temp)));
    l23_average_shape = temp/temp_max;
    shape_template=l23_average_shape;
    
    
    experiment_setup.prior_info.induced_intensity=get_spike_curves(experiment_setup.single_patch_path); 
    experiment_setup.prior_info.induced_intensity.minimum_stim_threshold = ...
       experiment_setup.prior_info.induced_intensity.current(max(find(experiment_setup.prior_info.induced_intensity.prob<0.2))+1);
    experiment_setup.prior_info.induced_intensity.stim_threshold=5;

    experiment_setup.neighbourhood_params=struct;
    
    experiment_setup.neighbourhood_params.z_bounds=[10 70]; % leave empty for no bounds
    experiment_setup.neighbourhood_params.x_bounds=[-150 150]; % leave empty for no bounds
    experiment_setup.neighbourhood_params.y_bounds=[-150 150]; % leave empty for no bounds
    experiment_setup.neighbourhood_params.number=10;
    experiment_setup.neighbourhood_params.height=15;
    experiment_setup.neighbourhood_params.buffer_height=40;  %NOTE: this count from the plane
    experiment_setup.neighbourhood_params.buffer_x=10;  % NOTE: These two from the boundariesedi
    experiment_setup.neighbourhood_params.buffer_y=10;
    
    experiment_setup.exp.z_width = 20;
    % experiment_setup.exp.z_depths = '10 30 50 70 90';% this should be a space delimited string
    experiment_setup.exp.arbitrary_z = 1;
    
    load([ experiment_setup.prior_root 'power-calibration.mat']);
    experiment_setup.exp.ratio_map = ratio_map;
    experiment_setup.exp.pockels_lut = pockels_lut;
    experiment_setup.exp.max_power_ref = max(experiment_setup.exp.pockels_lut(2,:));
    experiment_setup.exp.min_full_foe_power = experiment_setup.exp.max_power_ref/max(ratio_map(:));
    
    
    % experiment_setup.exp.min_iti = 55; % min iti for slidebook
    experiment_setup.exp.foe_bounds = [-148 148; -148 148; 0 400];
    
    % more experimental experiment_setup
    experiment_setup.exp.sim_locs = 0;
    experiment_setup.exp.stim_duration = .003; % length of laser on
    experiment_setup.run_parfor = 0;
    experiment_setup.exp.max_spike_freq = .5; % don't revisit cells on average sooner than this in Hz
    experiment_setup.exp.max_stim_freq = 1/(.050 + experiment_setup.exp.stim_duration + .001); % max frequency for system
    
    experiment_setup.phase_base_file = 1;
    if strcmp(experiment_setup.experiment_type,'experiment') || experiment_setup.sim.compute_phase_masks
        if experiment_setup.phase_base_file
            load('phase-mask-base-nocomp.mat');
        else % load from base ws
            %         disk_grid_phase = evalin('base','disk_grid_phase');
            disk_grid_key = evalin('base','disk_grid_key');
            %         fine_spots_grid_phase = evalin('base','fine_spots_grid_phase');
            fine_spots_grid_key = evalin('base','fine_spots_grid_key');
            
            % 2p image alignment params
            image_zero_order_coord = evalin('base','image_zero_order_coord');
            image_um_per_px = evalin('base','image_um_per_px');
            stack_um_per_slice = evalin('base','stack_um_per_slice');
            
        end
        %     experiment_setup.disk_grid_phase = disk_grid_phase;
        experiment_setup.disk_grid_key = disk_grid_key;
        %     experiment_setup.fine_spots_grid_phase = fine_spots_grid_phase;
        experiment_setup.fine_spots_grid_key = fine_spots_grid_key;
        
        experiment_setup.image_zero_order_coord = image_zero_order_coord;
        experiment_setup.image_um_per_px = image_um_per_px;
        experiment_setup.stack_um_per_slice = stack_um_per_slice;
        
    else
        %     experiment_setup.disk_grid_phase = [];
        experiment_setup.disk_grid_key = [];
        %     experiment_setup.fine_spots_grid_phase = [];
        experiment_setup.fine_spots_grid_key = [];
        
    end
    experiment_setup.exp.phase_mask_struct = 1;
    experiment_setup.exp.max_trials_per_sweep = 1000;
    experiment_setup.exp.first_stim_time = 1.0; % in sec
    experiment_setup.exp.filter_config = 'Femto Phasor';
    experiment_setup.exp.sweep_time_padding = 2.5; % in sec
    
    experiment_setup.exp.sim_response = 0;
%     experiment_setup.exp.run_online_detection = 1;
end

