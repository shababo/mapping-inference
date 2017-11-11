function [experiment_query_this_group] = design_undefined(this_neighbourhood,group_profile)
% Output:
% The experiment_query_this_group struct with fields:
%    experiment_query_this_group.trials(i_trial).location_IDs=this_trial_location_IDs;
%        experiment_query_this_group.trials(i_trial).cell_IDs=this_trial_cell_IDs;
%        experiment_query_this_group.trials(i_trial).power_levels=this_trial_power_levels;
%        experiment_query_this_group.trials(i_trial).locations=this_trial_locations;


group_ID=group_profile.group_ID;


%cells_this_group=find((arrayfun(@(x) strcmp(this_neighbourhood.neurons(1).group_ID,group_ID),this_neighbourhood.neurons(:))));
cells_this_group= find(get_group_inds(neighbourhood,group_ID));
number_cells_this_group=length(cells_this_group);
number_cells_all= length(this_neighbourhood.neurons);
loc_counts=zeros(number_cells_this_group,1);

% obtain posterior mean of gamma
% Write a function that grabs the last element in a specific field
i_batch=this_neighbourhood.batch_ID;
neurons=this_neighbourhood.neurons(cells_this_group);
properties={'PR_params'};
summary_stat={'mean'};
mean_gamma=grab_values_from_neurons(i_batch,neurons,properties,summary_stat);

this_neighbourhood.neurons(cells_this_group).PR_params(end).mean

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
        
        % calculate the firing probability
        firing_prob=cell(number_cells_this_group,1);
        
        for i_cell = 1:number_cells_this_group
            this_cell=cells_this_group(i_cell);
            candidate_grid=this_neighbourhood.neurons(this_cell).stim_locations.(group_ID);
            firing_prob{i_cell}=zeros(size(candidate_grid.effect,2),length(group_profile.design_func_params.trials_params.power_level),...
                size(candidate_grid.effect,1));
            for i_loc = 1:size(candidate_grid.effect,2)
                for k=1:length(group_profile.design_func_params.trials_params.power_level)
                    for j_cell = 1:size(candidate_grid.effect,2)
                        if candidate_grid.effect(j_cell,i_loc)>5e-2
                            stimulation_received=candidate_grid.effect(j_cell,i_loc)*group_profile.design_func_params.trials_params.power_level(k);
                            effective_stim= stimulation_received*gain_samples(:,j_cell);
                            stim_index=max(1,round(effective_stim*experiment_setup.prior_info.induced_intensity.stim_scale));
                            prob_collapsed=sum(experiment_setup.prior_info.induced_intensity.prob_trace_full(stim_index,:),2);
                            firing_prob{i_cell}(i_loc,k,j_cell)=mean(prob_collapsed);
                        end
                    end
                end
            end
        end
        loc_selected=zeros(number_cells_this_group,1);
        power_selected=zeros(number_cells_this_group,1);
        loc_to_cell_selected=1:number_cells_this_group;
        
        % Select the optimal locations based on firing_prob:
        for i_cell = 1:number_cells_this_group
            this_cell=cells_this_group(i_cell);
            firing_prob_temp=firing_prob{i_cell};
            firing_prob_temp(:,:,this_cell)=0;
            firing_prob_difference= firing_prob(:,:,this_cell)-max(firing_prob_temp,[],3);
            [max_value_loc,index_loc] = max(firing_prob_difference);
            % Pick the lowest power if the objectives are not too different from each
            % other
            weighted_max_value_loc = max_value_loc./log(group_profile.design_func_params.trials_params.power_level(k));
            weighted_max_value_loc( weighted_max_value_loc<0)=0;
            if max( weighted_max_value_loc)==0
                weighted_max_value_loc(:)=1;
            end
            index_I = ...
                randsample(1:length(weighted_max_value_loc),1,true,weighted_max_value_loc);
            %         [~,index_I]=max(weighted_max_value_loc);
            loc_selected(i_cell_idx)=index_loc(index_I);
            power_selected(i_cell_idx)=power_level(index_I);
        end
    case 'Nuclei'
        loc_selected=ones(number_cells_this_group,1);
        power_selected=zeros(number_cells_this_group,1);
    case 'Random'
        % do nothing here
    otherwise
        % throw a warning?
end

experiment_query_this_group=struct;
experiment_query_this_group.group_ID=group_ID;
experiment_query_this_group.neighbourhood_ID=this_neighbourhood.neighbourhood_ID;
experiment_query_this_group.instruction='Compute hologram';
experiment_query_this_group.trials=struct([]);
pockels_ratios = zeros(num_trials,n_spots_per_trial);

this_trial_location_IDs=zeros(1,group_profile.design_func_params.trials_params.spots_per_trial);
this_trial_cell_IDs=zeros(1,group_profile.design_func_params.trials_params.spots_per_trial);
this_trial_power_levels=zeros(1,group_profile.design_func_params.trials_params.spots_per_trial);
this_trial_locations=zeros(group_profile.design_func_params.trials_params.spots_per_trial,3);
for i_trial = 1:undefined_profile.design_func_params.trials_params.trials_per_batch
    prob_initial = probability_weights;
    prob_initial = prob_initial./(loc_counts+0.1);
    prob_initial = prob_initial/sum(prob_initial);
    pockels_ratio_refs(end+1) = 0;
    for i_spot = 1:undefined_profile.design_func_params.trials_params.spots_per_trial
        try_count = 0;
        loc_found = 0;
        if sum(prob_initial)>0.1
            prob_initial_thresh = prob_initial;
            thresh = median(prob_initial);
            prob_initial_thresh(prob_initial < thresh) = 0;
            while try_count < 10 && ~loc_found
                switch  group_profile.design_func_params.trials_params.stim_design
                    case 'Random'
                        temp_index = ...
                            randsample(1:number_cells_this_group,1,true,prob_initial_thresh);
                        this_cell=cells_this_group(temp_index);
                        grid_points_this_cell =size(this_neighbourhood.neurons(this_cell).stim_locations.(group_ID).grid,1);
                        temp_loc=randsample(1:grid_points_this_cell,1,true);
                    case 'Optimal'
                        temp_index = ...
                            randsample(1:number_cells_this_group,1,true,prob_initial_thresh);
                        this_cell=cells_this_group(temp_index);
                        
                        temp_loc =  loc_selected(temp_index);
                    case 'Nuclei'
                        temp_index = ...
                            randsample(1:number_cells_this_group,1,true,prob_initial_thresh);
                        this_cell=cells_this_group(temp_index);
                        
                        temp_loc =  loc_selected(temp_index);
                end
                
                
                %                     if use_power_map % NOT CODED FOR USE WITH REPLICATING WITHIN THIS FUNCTION!
                %                         this_loc = target_locations(temp_index,:);
                %                         ratio_this_loc = round(ratio_map(round(this_loc(1))+ceil(size(ratio_map,1)/2),...
                %                             round(this_loc(2))+ceil(size(ratio_map,2)/2))*10000);
                %                         total_ratio_tmp = pockels_ratio_refs(end) + ratio_this_loc/10000;
                %                         if ~(total_ratio_tmp > ratio_limit/power_selected(temp_loc)/i_spot) % this doesn't actually work for differnet powers per spot...
                %                             loc_found = 1;
                %                         end
                %                     else
                loc_found = 1;
                ratio_this_loc=0;
                total_ratio_tmp=0;
                %                     end
                try_count = try_count + 1;
            end
        end
        if ~loc_found
            this_trial_locations_ID(1,i_spot)=NaN;
            this_trial_power_levels(1,i_spot)=NaN;
        else
            loc_counts(temp_index)=loc_counts(temp_index)+1;
            %                 pockels_ratios(i_trial,i_spot) = ratio_this_loc;
            %                 pockels_ratio_refs(end) = total_ratio_tmp;
            this_trial_location_IDs(1,i_spot)=temp_loc;
            this_trial_cell_IDs(1,i_spot)=this_cell;
            switch  group_profile.design_func_params.trials_params.stim_design
                case 'Optimal'
                    this_trial_power_levels(1,i_spot)=power_selected(temp_index);
                    if rand(1) < mean_gamma(temp_index)
                        this_trial_powers(1,i_spot)=...
                            fire_stim_threshold./(pi_target(temp_index,loc_optimal(temp_index))...
                            *gain_samples(randsample(1:n_MC_samples,1),temp_index));
                        this_trial_powers(1,i_spot)=max(min(power_level),min(this_trial_powers(1,i_spot),max(power_level)));
                    end
                    
                otherwise
                    this_trial_power_levels(1,i_spot)=randsample(group_profile.design_func_params.trials_params.power_level,1,true);
            end
            this_trial_locations(i_spot,:)=    this_neighbourhood.neurons(this_cell).stim_locations.(group_ID).grid(temp_loc,:);
            prob_initial = prob_initial -  ...
                this_neighbourhood.neurons(this_cell).stim_locations.(group_ID).inner_normalized_products(cells_this_group,temp_loc)*prob_initial(temp_index);
            prob_initial = max(0,prob_initial);
            loc_found = 1;
        end
        
    end
    experiment_query_this_group.trials(i_trial).location_IDs=this_trial_location_IDs;
    experiment_query_this_group.trials(i_trial).cell_IDs=this_trial_cell_IDs;
    experiment_query_this_group.trials(i_trial).power_levels=this_trial_power_levels;
    experiment_query_this_group.trials(i_trial).locations=this_trial_locations;
    
end


end
