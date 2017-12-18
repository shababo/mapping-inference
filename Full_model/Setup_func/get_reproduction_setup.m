function reproduction_setup = get_reproduction_setup(varargin)
% General idea:
%  - Run each step of run_mapping_experiment (neighbourhood creations,
%  design trials, phase mask calibration, collect data (simulate data), 
%  event detection via OASIS, variational inference
%  - At each step, report discrepency between the reproduced results and
%  records, AND replace the reproduced results with records to continue the
%  process 


reproduction_setup=struct;
reproduction_setup.file_name='./Data/NewData/12_14_16_49/12_14_16_49_data.mat';

% FLAG FOR USING THE EXACT SAME SETTING
% NOTE: this should be relaxed to allow for using different hyperparameters
% etc for sensitivity analysis 
reproduction_setup.use_same_setting = true;

reproduction_setup.rep_params=struct;
reproduction_setup.rep_params.neighbourhoods = false; % Create neighbourhoods
reproduction_setup.rep_params.trials = false; % Design trials 
reproduction_setup.rep_params.phase_mask = false; % Calculate phase masks
reproduction_setup.rep_params.events = false;  % Simulate events
reproduction_setup.rep_params.VC = false;  % Simulate voltage clamp
reproduction_setup.rep_params.event_detection=false; % Run event detection
reproduction_setup.rep_params.inference=true; % Fit the model
reproduction_setup.rep_params.write_out=false; % whether to write out files

% FUNCTIONS FOR EVALUATING DISCREPENCY
reproduction_setup.rep_func=struct;
reproduction_setup.rep_func.neighbourhoods_comparison=[]; 
reproduction_setup.rep_func.designs_comparison=[]; 

% 
% 
% reproduction_setup.dimension=[-150 -150 0;150 150 30]; % there will be only two neighbourhoods 
% reproduction_setup.use_real_map=true;
% 
% reproduction_setup.real_map_name='6_3_s2c3_mpp_and_stim_data.mat'; % name of the cell map under 'experiment_setup.exp_root'
% reproduction_setup.number_of_cells=600;
% reproduction_setup.siblings=struct;
% reproduction_setup.siblings.number=20;
% reproduction_setup.siblings.distance=10;
% reproduction_setup.connection_params=struct;
% reproduction_setup.connection_params.proportion=0.2;
% reproduction_setup.cell_params=struct;
% reproduction_setup.cell_params.type='Normal';
% % Normal, Extreme gain, or Weak gamma.
% reproduction_setup.cell_params.gain_range=[0.5 0.9];
% reproduction_setup.cell_params.amp_mean = 2.5; %log norm
% reproduction_setup.cell_params.amp_std = .5; % log norm
% 
% reproduction_setup.do_instructions = 1;
% reproduction_setup.sim_vclamp = 0;
% reproduction_setup.response_generation=struct;
% 
% reproduction_setup.use_power_calib = 1;
% reproduction_setup.compute_phase_masks = 0;
% %       simulation_setup.response_generation.function=@xxx
%  
% reproduction_setup.visualize = 1;
% reproduction_setup.plotting_funcs = @plot_est_vs_truth;

