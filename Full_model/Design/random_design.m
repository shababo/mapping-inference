function [ trials_locations,  trials_powers] = random_design(...
    target_locations_selected,power_selected,...
    inner_normalized_products,single_spot_threshold,...
    gamma_estimates,prob_weight,...
    remaining_cell_list,n_spots_per_trial,K,n_replicates)

n_remaining_cell=length(remaining_cell_list);

if n_remaining_cell > single_spot_threshold
    n_unique_initial_per_plane = round(K*n_remaining_cell/n_spots_per_trial);
    trials_locations =zeros(n_unique_initial_per_plane*n_replicates,n_spots_per_trial);
    trials_powers =zeros(n_unique_initial_per_plane*n_replicates,n_spots_per_trial);
    n_stim_locations = size(target_locations_selected,1);
    this_trial_locations=zeros(1,n_spots_per_trial);
    this_trial_powers=zeros(1,n_spots_per_trial);
    % Set the probability to be inversely proportional to the
    % gamma_estimates, since we want to eliminate more disconnected cells
    probability_weights = (1-gamma_estimates)*prob_weight+(1-prob_weight); 
    
else
    
    n_unique_initial_per_plane = round(K*n_remaining_cell);
    this_trial_locations=zeros(1,n_spots_per_trial);
    this_trial_powers=zeros(1,n_spots_per_trial);
    this_trial_locations(1,2:end)=NaN;
    this_trial_powers(1,2:end)=NaN;
    n_spots_per_trial=1;
    probability_weights =ones(n_remaining_cell,1);
    
end
for i_trial =1:n_unique_initial_per_plane
    prob_initial = probability_weights;
    prob_initial=prob_initial/sum(prob_initial);
    for i_spot = 1:n_spots_per_trial
        if sum(prob_initial)>0.1
            temp_index = ...
                randsample(1:n_remaining_cell,1,true,prob_initial);
            temp_loc =  remaining_cell_list(temp_index);
            this_trial_locations(1,i_spot)=temp_loc;
            this_trial_powers(1,i_spot)=power_selected(temp_loc);
            prob_initial = ...
                prob_initial - inner_normalized_products(remaining_cell_list,temp_loc)*prob_initial(temp_index);
            prob_initial = max(0,prob_initial);
        else
            this_trial_locations(1,i_spot)=NaN;
            this_trial_powers(1,i_spot)=NaN;
        end
    end
    trials_locations(n_replicates*(i_trial-1)+(1:n_replicates),:)=ones(n_replicates,1)*this_trial_locations;
    trials_powers(n_replicates*(i_trial-1)+(1:n_replicates),:)=ones(n_replicates,1)*this_trial_powers;
end


end


