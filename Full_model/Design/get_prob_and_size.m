function [cells_probabilities, stimuli_size] = get_prob_and_size(...
    pi_target_selected,trials_locations,trials_powers,...
    stim_unique,prob_trace)
n_cell_this_plane =size(pi_target_selected,1);
stimuli_size=zeros(size(trials_locations,1),n_cell_this_plane);
for l = 1:size(trials_locations,1)
    for m = 1:size(trials_locations,2)
        if isnan(trials_locations(l,m))
        else
            stimuli_size(l,:) = stimuli_size(l,:)+...
                (pi_target_selected(:,trials_locations(l,m))*trials_powers(l,m))';
        end
    end
end
cells_probabilities = zeros([size(trials_locations,1) n_cell_this_plane]);
for   l = 1:size(trials_locations,1)
    for i_cell = 1:n_cell_this_plane
        cells_probabilities(l,i_cell)= stim_to_prob(stimuli_size(l,i_cell),stim_unique,prob_trace);
    end
end

end


