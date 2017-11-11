function [stimuli_size] = get_stim_size(group_ID,trials,this_neighbourhood)
    



number_of_trials=length(trials);
number_of_cells =length(this_neighbourhood.neurons);
stimuli_size=zeros(number_of_trials,number_of_cells);
number_of_spots = length(trials(1).location_IDs);
for i_trial = 1:number_of_trials
    for i_loc= 1:number_of_spots 
        if isnan(trials(i_trial).location_IDs(i_loc))
        else
            this_loc_ID=trials(i_trial).location_IDs(i_loc);
                this_loc_power = trials(i_trial).power_levels(i_loc);
                this_cell=trials(i_trial).cell_IDs(i_loc);
            stimuli_size(i_trial,:) = stimuli_size(i_trial,:)+...
                 (this_neighbourhood.neurons(this_cell).stim_locations.(group_ID).effect(:,this_loc_ID)*this_loc_power)';
             
        end
    end
end
end


