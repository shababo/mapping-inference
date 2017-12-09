function [experiment_query_this_group] = design_connected(this_neighbourhood,group_profile,experiment_setup)
% Output:
% The experiment_query_this_group struct with fields:
%    experiment_query_this_group.trials(i_trial).location_IDs=this_trial_location_IDs;
%        experiment_query_this_group.trials(i_trial).cell_IDs=this_trial_cell_IDs;
%        experiment_query_this_group.trials(i_trial).power_levels=this_trial_power_levels;
%        experiment_query_this_group.trials(i_trial).locations=this_trial_locations;


disp('designs connected')
group_ID=group_profile.group_ID;

power_levels=group_profile.design_func_params.trials_params.power_levels;
%i_cell_group_to_nhood=find((arrayfun(@(x) strcmp(this_neighbourhood.neurons(:).group_ID,group_ID),this_neighbourhood.neurons(:))));
i_cell_group_to_nhood= find(get_group_inds(this_neighbourhood,group_ID));
number_cells_this_group=length(i_cell_group_to_nhood);
number_cells_nhood= length(this_neighbourhood.neurons);
loc_counts=zeros(number_cells_this_group,1);

% obtain posterior mean of gamma
% Write a function that grabs the last element in a specific field
batch_ID=this_neighbourhood.batch_ID;
neurons_this_group=this_neighbourhood.neurons(i_cell_group_to_nhood);
properties={'PR_params'};summary_stat={'mean'};
temp_output=grab_values_from_neurons(batch_ID,neurons_this_group,properties,summary_stat);
mean_gamma=temp_output.PR_params.mean;

if  group_profile.design_func_params.trials_params.weighted_indicator
    probability_weights = 1-mean_gamma;
else
    probability_weights = ones(number_cells_this_group,1);
end

switch group_profile.design_func_params.trials_params.stim_design
    case 'Optimal'
        gain_samples=zeros(group_profile.inference_params.MCsamples_for_posterior,...
            number_cells_nhood);
        for i_cell_group = 1:number_cells_nhood
            alpha_gain=this_neighbourhood.neurons(i_cell_group).gain_params(end).alpha;
            beta_gain=this_neighbourhood.neurons(i_cell_group).gain_params(end).beta;
            temp=normrnd(alpha_gain,beta_gain,[group_profile.inference_params.MCsamples_for_posterior 1]);
            gain_samples(:,i_cell_group) = exp(temp)./(1+exp(temp))*...
                range(group_profile.inference_params.bounds.gain)+group_profile.inference_params.bounds.gain(1);
        end
        
        % calculate the firing probability
        firing_prob=cell(number_cells_this_group,1);
        
        for i_cell_group = 1:number_cells_this_group
            i_cell_nhood=i_cell_group_to_nhood(i_cell_group);
            candidate_grid=this_neighbourhood.neurons(i_cell_nhood).stim_locations.(group_ID);
            firing_prob{i_cell_group}=zeros(size(candidate_grid.effect,2),length(group_profile.design_func_params.trials_params.power_levels),...
                size(candidate_grid.effect,1));
            for i_loc = 1:size(candidate_grid.effect,2)
                for k=1:length(group_profile.design_func_params.trials_params.power_levels)
                    for j_cell = 1:size(candidate_grid.effect,1)
                        if candidate_grid.effect(j_cell,i_loc)>5e-2
                            stimulation_received=candidate_grid.effect(j_cell,i_loc)*group_profile.design_func_params.trials_params.power_levels(k);
                            effective_stim= stimulation_received*gain_samples(:,j_cell);
                            stim_index=max(1,round(effective_stim*experiment_setup.prior_info.induced_intensity.stim_scale));
                            prob_collapsed=sum(experiment_setup.prior_info.induced_intensity.intensity_grid(stim_index,:),2);
                            firing_prob{i_cell_group}(i_loc,k,j_cell)=mean(prob_collapsed);
                        end
                    end
                end
            end
        end
        loc_selected=zeros(number_cells_this_group,group_profile.design_func_params.trials_params.num_stim_sites);
        power_selected=zeros(number_cells_this_group,group_profile.design_func_params.trials_params.num_stim_sites);
        
        % Select the optimal locations based on firing_prob:
        for i_cell_group = 1:number_cells_this_group
            i_cell_nhood=i_cell_group_to_nhood(i_cell_group);
            candidate_grid=this_neighbourhood.neurons(i_cell_nhood).stim_locations.(group_ID).grid;
            
            firing_prob_temp=firing_prob{i_cell_group};
            firing_prob_temp(:,:,i_cell_nhood)=0;
            firing_prob_difference= firing_prob{i_cell_group}(:,:,i_cell_nhood)-max(firing_prob_temp,[],3);
            [max_value_loc,index_I] = max(firing_prob_difference');
            % Pick the lowest power if the objectives are not too different from each
            % other
            %weighted_max_value_loc = max_value_loc./log(group_profile.design_func_params.trials_params.power_levels);
            weighted_max_value_loc = max_value_loc;
            weighted_max_value_loc( weighted_max_value_loc<0)=0;
            if max( weighted_max_value_loc)==0
                weighted_max_value_loc(:)=1;
            end
            % for each cell find more places:
            still_available=weighted_max_value_loc;
            still_available(:)=1;
            %still_available
            for i_loc = 1:group_profile.design_func_params.trials_params.num_stim_sites
                %still_available
                %weighted_max_value_loc
                prob_available = weighted_max_value_loc.*still_available;
                
                if sum(prob_available) >0 % if there are still some options..
                if i_loc == 1
                    index_loc=1;
                else
                    index_loc = ...
                        randsample(1:length(weighted_max_value_loc),1,true,prob_available);
                end
                    % Update the availability:
%                     fprintf('%d',index_I);
                    na_index=find(sqrt(sum((candidate_grid - ones(size(candidate_grid,1),1)*candidate_grid(index_loc,:)).^2,2))<group_profile.design_func_params.trials_params.min_gap_stim);
                   still_available(na_index)=0;
                    
                else % Randomly pick one
                    index_loc = ...
                        randsample(1:length(weighted_max_value_loc),1,true,ones(length(weighted_max_value_loc),1));
                    
                end
                loc_selected(i_cell_group,i_loc)=index_loc;
                power_selected(i_cell_group,i_loc)=group_profile.design_func_params.trials_params.power_levels(index_I(index_loc));
                
            end
        end
        
    case 'Nuclei'
        
        loc_selected=ones(number_cells_this_group,1);
        power_selected=zeros(number_cells_this_group,1);
        
    case 'Random'
        
        loc_selected=zeros(number_cells_this_group,1);
        for i_cell_group = 1:number_cells_this_group
            i_cell_nhood=i_cell_group_to_nhood(i_cell_group);
            grid_points_this_cell = size(this_neighbourhood.neurons(i_cell_nhood).stim_locations.(group_ID).grid,1);
            loc_selected(i_cell_group)=randsample(1:grid_points_this_cell,1,true);
        end
        
        power_selected=zeros(number_cells_this_group,1);
        
    otherwise
        % throw a warning?
end

experiment_query_this_group=struct;
experiment_query_this_group.group_ID=group_ID;
experiment_query_this_group.neighbourhood_ID=this_neighbourhood.neighbourhood_ID;
experiment_query_this_group.instruction='Compute hologram';
experiment_query_this_group.trials=struct([]);
% pockels_ratios = zeros(num_trials,n_spots_per_trial);


%num_trials_per_cell = round(group_profile.design_func_params.trials_params.trials_per_batch*probability_weights/sum(probability_weights));
%num_trials_per_cell(num_trials_per_cell<group_profile.design_func_params.min_trials_per_cell)=group_profile.design_func_params.min_trials_per_cell;

num_trials_per_cell=ones(number_cells_this_group,1)*group_profile.design_func_params.trials_per_cell;
if sum(num_trials_per_cell)>group_profile.design_func_params.trials_params.trials_per_batch
    num_trials_per_cell= ceil( num_trials_per_cell*group_profile.design_func_params.trials_params.trials_per_batch/sum(num_trials_per_cell));
end
trial_counts=cumsum(num_trials_per_cell);
trial_counts=[0; trial_counts];

this_trial_location_IDs=zeros(1,group_profile.design_func_params.trials_params.spots_per_trial);
this_trial_cell_IDs=zeros(1,group_profile.design_func_params.trials_params.spots_per_trial);
this_trial_power_levels=zeros(1,group_profile.design_func_params.trials_params.spots_per_trial);
this_trial_locations=zeros(group_profile.design_func_params.trials_params.spots_per_trial,3);


for i_cell_group = 1:number_cells_this_group
    i_cell_nhood=i_cell_group_to_nhood(i_cell_group);
    for i_trial = (trial_counts(i_cell_group)+1):trial_counts(i_cell_group+1)
        this_trial_location_IDs=randsample(loc_selected(i_cell_group,:),1);
        i_cell_nhood=i_cell_group_to_nhood(i_cell_group);
        this_trial_cell_IDs= this_neighbourhood.neurons(i_cell_nhood).cell_ID;
        this_trial_power_levels=power_selected(i_cell_group);
        this_trial_locations=this_neighbourhood.neurons(i_cell_nhood).stim_locations.(group_ID).grid(this_trial_location_IDs,:);
        
        
        switch  group_profile.design_func_params.trials_params.stim_design
            case 'Optimal'
                this_trial_power_levels=power_selected(i_cell_group);
                if rand(1) < mean_gamma(i_cell_group)
                    this_trial_power_levels=...
                        experiment_setup.prior_info.induced_intensity.fire_stim_threshold./(neurons_this_group(i_cell_group).stim_locations.(group_ID).effect(i_cell_nhood,loc_selected(i_cell_group))...
                        *gain_samples(randsample(1:group_profile.inference_params.MCsamples_for_posterior,1),i_cell_nhood));
                    this_trial_power_levels=max(min(power_levels),min(this_trial_power_levels,max(power_levels)));
                end
                
            otherwise
                this_trial_power_levels=randsample(group_profile.design_func_params.trials_params.power_levels,1,true);
        end
        
        if strcmp(experiment_setup.experiment_type,'experiment') || experiment_setup.sim.use_power_calib
            high_power = 1;
            while high_power
                adj_pow_this_loc = round(experiment_setup.exp.ratio_map(round(this_trial_locations(1))+ceil(size(experiment_setup.exp.ratio_map,1)/2),...
                    round(this_trial_locations(2))+ceil(size(experiment_setup.exp.ratio_map,2)/2))*this_trial_power_levels);
                if adj_pow_this_loc > experiment_setup.exp.max_power_ref
                    this_trial_power_levels = this_trial_power_levels - 5;
                elseif this_trial_power_levels < experiment_setup.exp.min_full_foe_power
                    this_trial_power_levels = experiment_setup.exp.min_full_foe_power;
                    high_power = 0;
                else
                    high_power = 0;
                end
            end
        else
            adj_pow_this_loc = 0;
        end
        
        
        experiment_query_this_group.trials(i_trial).location_IDs=this_trial_location_IDs;
        experiment_query_this_group.trials(i_trial).cell_IDs=this_trial_cell_IDs;
        experiment_query_this_group.trials(i_trial).power_levels=this_trial_power_levels;
        experiment_query_this_group.trials(i_trial).locations=this_trial_locations;
        experiment_query_this_group.trials(i_trial).group_ID=group_ID;
        experiment_query_this_group.trials(i_trial).adj_power_per_spot = adj_pow_this_loc;
        
    end
    
end



