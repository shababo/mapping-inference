function [ trials_locations,  trials_powers] = random_design(...
    target_locations_selected,power_selected,...
    inner_normalized_products,single_spot_threshold,...
    gamma_estimates,prob_weight,...
    connected, loc_to_cell,...
    remaining_cell_list,n_spots_per_trial,K,n_replicates)

% connected: indicator whether we are designing trials for potentially
% connected cells 


n_remaining_cell=length(remaining_cell_list);
loc_counts=zeros(n_remaining_cell,1);
if n_remaining_cell < single_spot_threshold | n_spots_per_trial==1
    
    n_spots_per_trial=1;
    this_trial_locations=zeros(1,1);
    this_trial_powers=zeros(1,1);
    probability_weights =ones(n_remaining_cell,1);
    
    if connected
            
         for i_cell =1:n_remaining_cell % generate trials for each cell
            % Generate random locations on z-plane where the cell sits
            related_locations = find(loc_to_cell == i_cell);
            trials_temp=randsample(related_locations,K);
            
             trials_locations(n_replicates*K*(i_cell-1)+(1:(n_replicates*K)),:)=...
                 reshape([trials_temp; trials_temp],[n_replicates*K 1]);
            trials_powers(n_replicates*K*(i_cell-1)+(1:(n_replicates*K)),:)=...
                power_selected(trials_locations(n_replicates*K*(i_cell-1)+(1:(n_replicates*K)),:));
        end
        
    else % stim at the chosen locations
        
        for i_cell =1:n_remaining_cell % generate trials for each cell
            
            this_trial_locations=remaining_cell_list(i_cell);
            this_trial_powers=power_selected(this_trial_locations);
            trials_locations(n_replicates*K*(i_cell-1)+(1:(n_replicates*K)),:)=...
                ones(n_replicates*K,1)*this_trial_locations;
            trials_powers(n_replicates*K*(i_cell-1)+(1:(n_replicates*K)),:)=...
                ones(n_replicates*K,1)*this_trial_powers;
        end
        
    end

else
    n_unique_trials = round(K*n_remaining_cell/n_spots_per_trial);
    trials_locations =zeros(n_unique_trials *n_replicates,n_spots_per_trial);
    trials_powers =zeros(n_unique_trials*n_replicates,n_spots_per_trial);
    n_stim_locations = size(target_locations_selected,1);
    this_trial_locations=zeros(1,n_spots_per_trial);
    this_trial_powers=zeros(1,n_spots_per_trial);
    % Set the probability to be inversely proportional to the
    % gamma_estimates, since we want to eliminate more disconnected cells
    probability_weights = (1-gamma_estimates)*prob_weight+(1-prob_weight);
    
    for i_trial =1:n_unique_trials
        prob_initial = probability_weights;
        prob_initial=prob_initial./(loc_counts+0.1);
        prob_initial=prob_initial/sum(prob_initial);
        for i_spot = 1:n_spots_per_trial
            if sum(prob_initial)>0.1
                temp_index = ...
                    randsample(1:n_remaining_cell,1,true,prob_initial);
                loc_counts(temp_index)=loc_counts(temp_index)+1;
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

end


