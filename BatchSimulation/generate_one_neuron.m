function [one_neuron] = generate_one_neuron(simulation_params,prior_info)
one_neuron=struct;
    one_neuron.truth=struct;
%     gamma_truth = (rand([number_of_cells 1])<simulation_params.connection_params.proportion).*(0.5+0.3*rand([number_of_cells 1]));
        
%     one_neuron.truth.location=cell_locations(i_cell,:);
    shift_tmp=[0 0 0];
    one_neuron.truth.shift=shift_tmp;
    one_neuron.truth.optical_gain=0.02+(rand([1 1])-0.5)*0.01;
%     one_neuron.truth.PR=gamma_truth(i_cell);
    
    one_neuron.truth.delay_mean=(rand(1)-0.5)*40+20;
    one_neuron.truth.delay_var=(rand(1)-0.5)*20+15;
    one_neuron.truth.shape=[];
    one_neuron.fluorescence= []; % need to generate some fluorescence level
    
%     neurons(i_cell).cell_ID=i_cell;
    
%     neurons(i_cell).location=neurons(i_cell).truth.location+neurons(i_cell).truth.shift;
% Draw shapes:
GP_params=prior_info.prior_parameters.GP_params;
tmp= draw_3D_GP(simulation_params.mesh_grid,1,GP_params);
one_neuron.truth.shape= tmp.full.samples(:,1); 
