function [stimuli_size] = get_stim_size(group_ID,trials,this_neighbourhood)
    



number_of_trials=length(trials);
number_of_cells =length(this_neighbourhood.neurons);
number_of_spots = length(trials(1).location_IDs);

stimuli_size=zeros(number_of_trials,number_of_cells);

for l = 1:number_of_trials
    for m = 1:number_of_spots
        this_loc_ID=experiment_query_this_group.trials(l).location_IDs(m);
        if isnan(this_loc_ID)
        else
            this_loc_power = experiment_query_this_group.trials(l).power_levels(m);
            cell_ID=experiment_query_this_group.trials(l).cell_IDs(m);
            stimuli_size(l,:) = stimuli_size(l,:)+...
                (this_neighbourhood.neurons(cell_ID_list == cell_ID).stim_locations.(group_ID).effect(:,this_loc_ID)*this_loc_power)';
        end
    end
end


