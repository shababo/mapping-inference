function [ trials_locations,  trials_powers, target_locations_key, pockels_ratio_refs, pockels_ratios] = random_design(...
    target_locations_selected,power_selected,...
    inner_normalized_products,single_spot_threshold,...
    gamma_estimates,prob_weight,...
    connected, loc_to_cell,...
    remaining_cell_list,n_spots_per_trial,K,n_replicates,...
    varargin)


% connected: indicator whether we are designing trials for potentially
% connected cells 

if ~isempty(varargin) && ~isempty(varargin{1})
    use_power_map = varargin{1};
else
    use_power_map = 0;
end

if use_power_map
    ratio_map = varargin{2};
    ratio_limit = varargin{3};
end

if length(varargin) > 4 && ~isempty(varargin{4})
    do_replicates = varargin{4};
else
    do_replicates = 1;
end

if ~do_replicates
    n_replicates = 1;
%     K = 1;
end

pockels_ratio_refs = [];

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
            related_locations = find(loc_to_cell == remaining_cell_list(i_cell) );
            trials_temp=randsample(related_locations,K,true);
            
            trials_locations(n_replicates*K*(i_cell-1)+(1:(n_replicates*K)),:)=...
                 reshape(repmat(trials_temp,1,n_replicates),[n_replicates*K 1]);
            trials_powers(n_replicates*K*(i_cell-1)+(1:(n_replicates*K)),:)=...
                power_selected(trials_locations(n_replicates*K*(i_cell-1)+(1:(n_replicates*K)),:));
            
            if use_power_map 
                for i = 1:length(trials_temp)
                    for j = 1:n_replicates
                        this_loc = target_locations_selected(trials_temp(i),:);
                        pockels_ratio_refs(end + 1) = round(ratio_map(round(this_loc(1))+ceil(size(ratio_map,1)/2),...
                                 round(this_loc(2))+ceil(size(ratio_map,2)/2))*10000)/10000;
                    end
                end
            end
        end
        
    else % stim at the chosen locations
        
        for i_cell =1:n_remaining_cell % generate trials for each cell
            
            this_trial_locations=remaining_cell_list(i_cell);
            this_trial_powers=power_selected(this_trial_locations);
            trials_locations(n_replicates*K*(i_cell-1)+(1:(n_replicates*K)),:)=...
                ones(n_replicates*K,1)*this_trial_locations;
            trials_powers(n_replicates*K*(i_cell-1)+(1:(n_replicates*K)),:)=...
                ones(n_replicates*K,1)*this_trial_powers;
            if use_power_map 
                for i = 1:length(this_trial_locations)
                    for j = 1:n_replicates*K
                        this_loc = target_locations_selected(this_trial_locations(i),:);
                        pockels_ratio_refs(end + 1) = round(ratio_map(round(this_loc(1))+ceil(size(ratio_map,1)/2),...
                                 round(this_loc(2))+ceil(size(ratio_map,2)/2))*10000)/10000;
                    end
                end
            end
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
    pockels_ratios = zeros(n_unique_trials,n_spots_per_trial);
    for i_trial = 1:n_unique_trials
        prob_initial = probability_weights;
        prob_initial = prob_initial./(loc_counts+0.1);
        prob_initial = prob_initial/sum(prob_initial);
        pockels_ratio_refs(end+1) = 0
        for i_spot = 1:n_spots_per_trial
            
           if use_power_map
               power_test = pockels_ratio_refs(i_trial) < ratio_limit;
           else
               power_test = 1;
           end
           
            if sum(prob_initial)>0.1 && power_test
                temp_index = ...
                    randsample(1:n_remaining_cell,1,true,prob_initial);
                if use_power_map % NOT CODED FOR USE WITH REPLICATING WITHIN THIS FUNCTION!
                    this_loc = target_locations_selected(temp_index,:);
                    pockels_ratios(i_trial,i_spot) = round(ratio_map(round(this_loc(1))+ceil(size(ratio_map,1)/2),...
                                                                 round(this_loc(2))+ceil(size(ratio_map,2)/2))*10000);
                    pockels_ratio_refs(i_trial) = pockels_ratio_refs(i_trial) + pockels_ratios(i_trial,i_spot)/10000;
                    if pockels_ratio_refs(i_trial) > ratio_limit
                        this_trial_locations(1,i_spot)=NaN;
                        this_trial_powers(1,i_spot)=NaN;
                        continue
                    end
                end
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

target_locations_key = zeros(size(trials_locations,1),3,n_spots_per_trial);
for i = 1:size(trials_locations,1)
    for k = 1:n_spots_per_trial
        target_ind = trials_locations(i,k);
        target_locations_key(i,:,k) = target_locations_selected(target_ind,:);
    end
end

end


