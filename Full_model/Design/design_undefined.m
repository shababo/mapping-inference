function [experiment_query_this_group] = design_undefined(this_neighbourhood,group_profile)

group_type_ID=group_profile.group_type_ID;
cells_this_group=index([this_neighbourhoods.neurons(:).group_type_ID]==group_type_ID);
number_cells_this_group=length(cells_this_group);
number_cells_all= length(this_neighbourhoods.neurons);
loc_counts=zeros(number_cells_this_group,1);

% obtain posterior mean of gamma
% Write a function that grabs the last element in a specific field
mean_gamma=grab_recent_value(this_neighbourhoods.neurons(cells_this_group).PR_params(end).mean);
if  group_profile.design_func_params.trials_params.weighted_indicator
   probability_weights = 1-mean_gamma;
else
    probability_weights =ones(number_cells_this_group,1);
end
   
switch group_profile.design_func_params.trials_params.stim_design
    case 'Optimal'
        
        gain_samples=zeros(undefined_profile.inference_params.MCsamples_for_posterior,...
            number_cells_all);
        for i_cell = 1:number_cells_all
            alpha_gain=this_neighbourhood.neurons(i_cell).gain_params(end).alpha;
            beta_gain=this_neighbourhood.neurons(i_cell).gain_params(end).beta;
            temp=normrnd(alpha_gain,beta_gain,[undefined_profile.inference_params.MCsamples_for_posterior 1]);
            gain_samples(:,i_cell) = exp(temp)./(1+exp(temp))*...
                range(group_profile.inference_params.bounds.gain)+group_profile.inference_params.bounds.gain(1);
        end
        % need a function that creates the pi_target:
        pi_target=
        % calculate the firing probability
        firing_prob=zeros(size(pi_target,2),length(group_profile.design_func_params.trials_params.power_level),...
            size(pi_target,1));
        for i_loc = 1:size(pi_target,2)
            for k=1:length(power_level)
                for i_cell = 1:size(pi_target,1)
                    if pi_target(i_cell,i_loc)>5e-2
                        stimulation_received=pi_target(i_cell,i_loc)*power_level(k);
                        effective_stim= stimulation_received*gain_samples(:,i_cell);
                        stim_index=max(1,round(effective_stim*stim_scale));
                        prob_collapsed=sum(prob_trace_full(stim_index,:),2);
                        firing_prob(i_loc,k,i_cell)=mean(prob_collapsed);
                    end
                end
            end
        end
        loc_optimal=zeros(n_remaining_cell,1);
        power_optimal=zeros(n_remaining_cell,1);
        loc_to_cell_optimal=1:n_remaining_cell;
        
        % Select the optimal locations based on firing_prob:
        for i_cell_idx = 1:n_remaining_cell
            i_cell=remaining_cell_list(i_cell_idx);
            firing_prob_temp=firing_prob;
            firing_prob_temp(:,:,i_cell)=0;
            firing_prob_difference= firing_prob(:,:,i_cell)-max(firing_prob_temp,[],3);
            [max_value_loc,index_loc] = max(firing_prob_difference);
            % Pick the lowest power if the objectives are not too different from each
            % other
            weighted_max_value_loc = max_value_loc./log(power_level);
            weighted_max_value_loc( weighted_max_value_loc<0)=0;
            if max( weighted_max_value_loc)==0
                weighted_max_value_loc(:)=1;
            end
            index_I = ...
                randsample(1:length(weighted_max_value_loc),1,true,weighted_max_value_loc);
            %         [~,index_I]=max(weighted_max_value_loc);
            loc_optimal(i_cell_idx)=index_loc(index_I);
            power_optimal(i_cell_idx)=power_level(index_I);
        end
    case 'Nuclei'
    case 'Random'
    otherwise
        %a warning?
end

    if design_type_multi == 1 % optimal stimulation
        % New design for multi-spot stimulation:
        % Draw samples from the current posterior distribution of gains
       
    end
    
    % Design trials using the optimal locations & power
    
    trials_locations =zeros(num_trials,n_spots_per_trial);
    trials_powers =zeros(num_trials,n_spots_per_trial);
    n_stim_locations = size(target_locations,1);
    this_trial_locations=zeros(1,n_spots_per_trial);
    this_trial_powers=zeros(1,n_spots_per_trial);
    
    % Set the probability to be inversely proportional to the
    % gamma_estimates, since we want to eliminate more disconnected cells
    pockels_ratios = zeros(num_trials,n_spots_per_trial);
    for i_trial = 1:num_trials
        prob_initial = probability_weights;
        prob_initial = prob_initial./(loc_counts+0.1);
        prob_initial = prob_initial/sum(prob_initial);
        pockels_ratio_refs(end+1) = 0;
        for i_spot = 1:n_spots_per_trial
            
            %            if use_power_map
            %                power_test = pockels_ratio_refs(end) < ratio_limit/power_selected(1)/i_spot; % hack here since all same power...
            %            else
            %                power_test = 1;
            %            end
            try_count = 0;
            loc_found = 0;
            if sum(prob_initial)>0.1
                prob_initial_thresh = prob_initial;
                thresh = median(prob_initial);
                prob_initial_thresh(prob_initial < thresh) = 0;
                while try_count < 10 && ~loc_found
                    switch design_type_multi
                        case 0 % random stimulation
                        temp_index = ...
                            randsample(1:n_remaining_cell,1,true,prob_initial_thresh);
                        related_locations = find(loc_to_cell == cell_list_map(remaining_cell_list(temp_index)));
                        temp_loc=randsample(related_locations,1,true);
                        
                        case 1 % optimized stimulation
                        temp_index = ...
                            randsample(1:n_remaining_cell,1,true,prob_initial_thresh);
                        temp_loc =  loc_optimal(temp_index);
                        case 2 % nuclei only
                        temp_index = ...
                            randsample(1:n_remaining_cell,1,true,prob_initial_thresh);
                        related_locations = find(loc_to_cell == cell_list_map(remaining_cell_list(temp_index)));
                        temp_loc=related_locations(1);
                    end
                    
                    
                    if use_power_map % NOT CODED FOR USE WITH REPLICATING WITHIN THIS FUNCTION!
                        this_loc = target_locations(temp_index,:);
                        ratio_this_loc = round(ratio_map(round(this_loc(1))+ceil(size(ratio_map,1)/2),...
                            round(this_loc(2))+ceil(size(ratio_map,2)/2))*10000);
                        total_ratio_tmp = pockels_ratio_refs(end) + ratio_this_loc/10000;
                        if ~(total_ratio_tmp > ratio_limit/power_selected(temp_loc)/i_spot) % this doesn't actually work for differnet powers per spot...
                            loc_found = 1;
                        end
                    else
                        loc_found = 1;
                        ratio_this_loc=0; % so that the code can run if use_power_map =0;
                        total_ratio_tmp=0;
                    end
                    try_count = try_count + 1;
                end
            end
            if ~loc_found
                this_trial_locations(1,i_spot)=NaN;
                this_trial_powers(1,i_spot)=NaN;
            else
                loc_counts(temp_index)=loc_counts(temp_index)+1;
                pockels_ratios(i_trial,i_spot) = ratio_this_loc;
                pockels_ratio_refs(end) = total_ratio_tmp;
                this_trial_locations(1,i_spot)=temp_loc;
                
                if design_type_single == 1 % optimal stimulation
                    
                    this_trial_powers(1,i_spot)=power_optimal(temp_index);
                
                
                % with probability proportional to its gamma, jitter the
                % power:
                if rand(1) < mean_gamma(temp_index)
                    this_trial_powers(1,i_spot)=...
                        fire_stim_threshold./(pi_target(temp_index,loc_optimal(temp_index))...
                        *gain_samples(randsample(1:n_MC_samples,1),temp_index));
                    this_trial_powers(1,i_spot)=max(min(power_level),min(this_trial_powers(1,i_spot),max(power_level)));
                end
                    
                else % otherwise
                    this_trial_powers(1,i_spot)=randsample(power_level,1,true);
                
                end
                
                prob_initial = ...
                    prob_initial - inner_normalized_products(remaining_cell_list,temp_loc)*prob_initial(temp_index);
                prob_initial = max(0,prob_initial);
                loc_found = 1;
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
        if isnan(target_ind)
            target_locations_key(i,:,k) = NaN;
        else
            target_locations_key(i,:,k) = target_locations(target_ind,:);
        end
    end
end

end