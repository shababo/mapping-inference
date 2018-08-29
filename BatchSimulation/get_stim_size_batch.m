function stimuli_size = get_stim_size_batch(neurons,trials,simulation_params)
number_of_trials = length(trials);
number_of_cells =length(neurons);
stimuli_size=zeros(number_of_trials,number_of_cells);
for l=1:number_of_trials
    this_loc=trials(l).locations; % 1 spot per trial
    this_power=trials(l).power_levels; % 1 spot per trial
    for i_cell=1:number_of_cells
        rel_loc =  this_loc - neurons(i_cell).truth.location;
        this_size = griddata(simulation_params.mesh_grid(:,1),simulation_params.mesh_grid(:,2),simulation_params.mesh_grid(:,3),...
            neurons(i_cell).truth.shape,rel_loc(1),rel_loc(2),rel_loc(3));
        if isnan(this_size)
            this_size = 0;
        end
        stimuli_size(l,i_cell) = this_power*this_size;
    end
end