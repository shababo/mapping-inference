function simulation_setup = get_simulation_setup(varargin)
% to-dos:
% 1) allow for generating siblings in xy-plain or z-axis
% 2) revise the plot function plot_est_vs_truth
simulation_setup=struct;


simulation_setup.dimension=[-50 -50 0;50 50 30]; % there will be only two neighbourhoods 
simulation_setup.use_real_map=true;
simulation_setup.background_rate= 1e-5;
simulation_setup.real_map_name='6_3_s2c3_mpp_and_stim_data.mat'; % name of the cell map under 'experiment_setup.exp_root'
simulation_setup.number_of_cells=10;
simulation_setup.siblings=struct;
simulation_setup.siblings.number=1;
simulation_setup.siblings.distance=15; 
simulation_setup.connection_params=struct;
simulation_setup.connection_params.proportion=0.2;
simulation_setup.cell_params=struct;
simulation_setup.cell_params.type='Normal';
% Normal, Extreme gain, or Weak gamma.
% simulation_setup.cell_params.gain_range=[0.5 0.9];
% simulation_setup.cell_params.amp_mean = 2.5; %log norm
% simulation_setup.cell_params.amp_std = .5; % log norm

simulation_setup.do_instructions = 0;
simulation_setup.sim_vclamp = 0;
simulation_setup.response_generation=struct;

simulation_setup.use_power_calib = 1;
simulation_setup.compute_phase_masks = 0;
%       simulation_setup.response_generation.function=@xxx
 
simulation_setup.visualize = 1;
simulation_setup.plotting_funcs = @plot_est_vs_truth;
simulation_setup.grid.z = -80:5:80;
simulation_setup.grid.x = -26:4:26;
simulation_setup.grid.y = -26:4:26;
[tmp_x, tmp_y, tmp_z]= meshgrid(simulation_setup.grid.x,simulation_setup.grid.y,simulation_setup.grid.z);
all_dim = length(simulation_setup.grid.z)*length(simulation_setup.grid.y)*length(simulation_setup.grid.x);
simulation_setup.mesh_grid = [reshape(tmp_x, [all_dim 1]) reshape(tmp_y, [all_dim 1]) reshape(tmp_z, [all_dim 1])];


