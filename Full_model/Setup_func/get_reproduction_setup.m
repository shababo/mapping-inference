function reproduction_setup = get_reproduction_setup(varargin)


reproduction_setup=struct;
reproduction_setup.file_name='./Data/NewData/12_14_16_49/12_14_16_49_data.mat';

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

