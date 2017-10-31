function [stimuli_size] = get_stim_size(...
    pi_target,trials_locations,trials_powers)
n_cell_this_plane =size(pi_target,1);
stimuli_size=zeros(size(trials_locations,1),n_cell_this_plane);
for l = 1:size(trials_locations,1)
    for m = 1:size(trials_locations,2)
        if isnan(trials_locations(l,m))
        else
            stimuli_size(l,:) = stimuli_size(l,:)+...
                (pi_target(:,trials_locations(l,m))*trials_powers(l,m))';
        end
    end
end
end


