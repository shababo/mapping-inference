
function simulation_setup = get_simulation_setup(varargin)


simulation_setup=struct;
simulation_setup.dimension=[-150 150;-150 150;0 150];
simulation_setup.number_of_cells=500;
simulation_setup.siblings=struct;
simulation_setup.siblings.number=500;
simulation_setup.siblings.distance=10;
simulation_setup.connection_params=struct;
simulation_setup.connection_params.proportion=0.2;
simulation_setup.cell_params=struct;
simulation_setup.cell_params.type='Normal';
% Normal, Extreme gain, or Weak gamma.
simulation_setup.cell_params.gain_range=[0.5 0.9];


simulation_setup.sim_vclamp = 0;
        