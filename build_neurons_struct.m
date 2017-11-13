function neurons = build_neurons_struct(nuclear_locs,fluor_vals,experiment_setup)

number_of_cells = size(nuclear_locs,1);


neurons=struct([]);
for i_cell = 1:number_of_cells 
    neurons(i_cell).fluorescence= fluor_vals(i_cell); % need to generate some fluorescence level 
    neurons(i_cell).V_reset= -1e4;
    neurons(i_cell).V_thresh=15;
    neurons(i_cell).membrane_resistance=0.02;
    neurons(i_cell).location=nuclear_locs(i_cell,1:3);
    neurons(i_cell).optical_gain=[];
    neurons(i_cell).PR=[];
    neurons(i_cell).truth=struct;
        neurons(i_cell).truth.V_reset= -1e4;
        neurons(i_cell).truth.V_thresh=15;
        neurons(i_cell).truth.membrane_resistance=0.02;
        neurons(i_cell).truth.location=nuclear_locs(i_cell,:);
        neurons(i_cell).truth.optical_gain=NaN;
        neurons(i_cell).truth.PR=NaN;
end