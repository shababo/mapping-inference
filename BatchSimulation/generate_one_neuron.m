function [one_neuron] = generate_one_neuron(simulation_params)
one_neuron=struct;
    one_neuron.truth=struct;
%     gamma_truth = (rand([number_of_cells 1])<simulation_params.connection_params.proportion).*(0.5+0.3*rand([number_of_cells 1]));
        
%     one_neuron.truth.location=cell_locations(i_cell,:);
    shift_tmp=[0 0 0];
    one_neuron.truth.shift=shift_tmp;
    one_neuron.truth.optical_gain=0.02+(rand([1 1])-0.5)*0.01;
%     one_neuron.truth.PR=gamma_truth(i_cell);
    
    one_neuron.truth.delay_mean=(rand(1)-0.5)*5+40;
    one_neuron.truth.delay_var=(rand(1)-0.5)*20+18;
    one_neuron.truth.shape=[];
        one_neuron.truth.PR=0.8;
        
        one_neuron.truth.background=simulation_params.background_rate;
     one_neuron.truth.location=[0 0 0];
    one_neuron.fluorescence= []; % need to generate some fluorescence level
 one_neuron.location=[0 0 0];
 one_neuron.cell_ID=1;
 
    