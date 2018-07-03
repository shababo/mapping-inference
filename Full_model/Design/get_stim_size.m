function [stimuli_size] = get_stim_size(trials,this_neighbourhood,experiment_setup)
%% Debug:
% this_neighbourhood = neighbourhoods(1);
% trials = experiment_query.undefined.trials;
%% Dimensions 
number_of_trials=length(trials);
number_of_cells =length(this_neighbourhood.neurons);
number_of_spots = size(trials(1).locations,2);
simulation_params=experiment_setup.sim;
%%
stimuli_size=zeros(number_of_trials,number_of_cells);
cell_ID_list=[this_neighbourhood.neurons(:).cell_ID];

for l = 1:number_of_trials
    for m = 1:number_of_spots
        this_loc_ID=trials(l).location_IDs(m);
        this_loc=trials(l).locations(m,:);
        if isnan(this_loc_ID)
        else
            this_loc_power =trials(l).power_levels(m);
            cell_ID=trials(l).cell_IDs(m);
            i_cell = find(cell_ID_list == cell_ID);
            rel_loc = this_loc - this_neighbourhood.neurons(i_cell).truth.location;
            if this_neighbourhood.neurons(i_cell).truth.PR>0
            stimuli_size(l,:) = stimuli_size(l,:)+...
                this_loc_power*griddata(simulation_params.mesh_grid(:,1),simulation_params.mesh_grid(:,2),simulation_params.mesh_grid(:,3),...
                this_neighbourhood.neurons(i_cell).truth.shape,rel_loc(1),rel_loc(2),rel_loc(3));
            end
        end
    end
end


