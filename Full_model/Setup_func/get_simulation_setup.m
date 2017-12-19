function simulation_setup = get_simulation_setup(varargin)


simulation_setup=struct;


simulation_setup.dimension=[-160 -160 0;160 160 80]; % there will be only two neighbourhoods 
simulation_setup.use_real_map=true;

simulation_setup.real_map_name='6_3_s2c3_mpp_and_stim_data.mat'; % name of the cell map under 'experiment_setup.exp_root'
simulation_setup.number_of_cells=450;
simulation_setup.siblings=struct;
simulation_setup.siblings.number=20;
simulation_setup.siblings.distance=10;
simulation_setup.connection_params=struct;
simulation_setup.connection_params.proportion=0.2;
simulation_setup.cell_params=struct;
simulation_setup.cell_params.type='Normal';
% Normal, Extreme gain, or Weak gamma.
simulation_setup.cell_params.gain_range=[0.5 0.9];
simulation_setup.cell_params.amp_mean = 2.5; %log norm
simulation_setup.cell_params.amp_std = .5; % log norm

simulation_setup.do_instructions = 0;
simulation_setup.sim_vclamp = 0;
simulation_setup.response_generation=struct;

simulation_setup.use_power_calib = 1;
simulation_setup.compute_phase_masks = 0;
%       simulation_setup.response_generation.function=@xxx
 
simulation_setup.visualize = 1;
simulation_setup.plotting_funcs = @plot_est_vs_truth;

