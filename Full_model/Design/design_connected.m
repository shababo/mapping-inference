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
%cells_this_group=find((arrayfun(@(x) strcmp(this_neighbourhood.neurons(:).group_ID,group_ID),this_neighbourhood.neurons(:))));
cells_this_group= find(get_group_inds(this_neighbourhood,group_ID));
number_cells_this_group=length(cells_this_group);
number_cells_all= length(this_neighbourhood.neurons);
loc_counts=zeros(number_cells_this_group,1);

% obtain posterior mean of gamma
% Write a function that grabs the last element in a specific field
i_batch=this_neighbourhood.batch_ID;
neurons=this_neighbourhood.neurons(cells_this_group);
properties={'PR_params'};summary_stat={'mean'};
temp_output=grab_values_from_neurons(i_batch,neurons,properties,summary_stat);
mean_gamma=temp_output.PR_params.mean;

if  group_profile.design_func_params.trials_params.weighted_indicator
    probability_weights = 1-mean_gamma;
else
    probability_weights =ones(number_cells_this_group,1);
end

switch group_profile.design_func_params.trials_params.stim_design
    case 'Optimal'
        gain_samples=zeros(group_profile.inference_params.MCsamples_for_posterior,...
            number_cells_all);
        for i_cell = 1:number_cells_all
            alpha_gain=this_neighbourhood.neurons(i_cell).gain_params(end).alpha;
            beta_gain=this_neighbourhood.neurons(i_cell).gain_params(end).beta;
            temp=normrnd(alpha_gain,beta_gain,[group_profile.inference_params.MCsamples_for_posterior 1]);
            gain_samples(:,i_cell) = exp(temp)./(1+exp(temp))*...
                range(group_profile.inference_params.bounds.gain)+group_profile.inference_params.bounds.gain(1);
        end
        
        % calculate the firing probability
        firing_prob=cell(number_cells_this_group,1);
        
        for i_cell = 1:number_cells_this_group
            this_cell=cells_this_group(i_cell);
            candidate_grid=this_neighbourhood.neurons(this_cell).stim_locations.(group_ID);
            firing_prob{i_cell}=zeros(size(candidate_grid.effect,2),length(group_profile.design_func_params.trials_params.power_levels),...
                size(candidate_grid.effect,1));
            for i_loc = 1:size(candidate_grid.effect,2)
                for k=1:length(group_profile.design_func_params.trials_params.power_levels)
                    for j_cell = 1:size(candidate_grid.effect,1)
                        if candidate_grid.effect(j_cell,i_loc)>5e-2
                            stimulation_received=candidate_grid.effect(j_cell,i_loc)*group_profile.design_func_params.trials_params.power_levels(k);
                            effective_stim= stimulation_received*gain_samples(:,j_cell);
                            stim_index=max(1,round(effective_stim*experiment_setup.prior_info.induced_intensity.stim_scale));
                            prob_collapsed=sum(experiment_setup.prior_info.induced_intensity.intensity_grid(stim_index,:),2);
                            firing_prob{i_cell}(i_loc,k,j_cell)=mean(prob_collapsed);
                        end
                    end
                end
            end
        end
        loc_selected=zeros(number_cells_this_group,1);
        power_selected=zeros(number_cells_this_group,1);
        
        % Select the optimal locations based on firing_prob:
        for i_cell = 1:number_cells_this_group
            this_cell=cells_this_group(i_cell);
            firing_prob_temp=firing_prob{i_cell};
            firing_prob_temp(:,:,this_cell)=0;
            firing_prob_difference= firing_prob{i_cell}(:,:,this_cell)-max(firing_prob_temp,[],3);
            [max_value_loc,index_loc] = max(firing_prob_difference);
            % Pick the lowest power if the objectives are not too different from each
            % other
            weighted_max_value_loc = max_value_loc./log(group_profile.design_func_params.trials_params.power_levels(k));
            weighted_max_value_loc( weighted_max_value_loc<0)=0;
            if max( weighted_max_value_loc)==0
                weighted_max_value_loc(:)=1;
            end
            index_I = ...
                randsample(1:length(weighted_max_value_loc),1,true,weighted_max_value_loc);
            %         [~,index_I]=max(weighted_max_value_loc);
            loc_selected(i_cell)=index_loc(index_I);
            power_selected(i_cell)=group_profile.design_func_params.trials_params.power_levels(index_I);
        end
          
    case 'Nuclei'
         
        loc_selected=ones(number_cells_this_group,1);
        power_selected=zeros(number_cells_this_group,1);
    case 'Random'
        loc_selected=zeros(number_cells_this_group,1); 
        for i_cell = 1:number_cells_this_group
            grid_points_this_cell =size(this_neighbourhood.neurons(this_cell).stim_locations.(group_ID).grid,1);
            loc_selected(i_cell)=randsample(1:grid_points_this_cell,1,true);
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


num_trials_per_cell = round(group_profile.design_func_params.trials_params.trials_per_batch*probability_weights/sum(probability_weights));
num_trials_per_cell(num_trials_per_cell<group_profile.design_func_params.min_trials_per_cell)=group_profile.design_func_params.min_trials_per_cell;
trial_counts=cumsum(num_trials_per_cell);
trial_counts=[0; trial_counts];

this_trial_location_IDs=zeros(1,group_profile.design_func_params.trials_params.spots_per_trial);
this_trial_cell_IDs=zeros(1,group_profile.design_func_params.trials_params.spots_per_trial);
this_trial_power_levels=zeros(1,group_profile.design_func_params.trials_params.spots_per_trial);
this_trial_locations=zeros(group_profile.design_func_params.trials_params.spots_per_trial,3);


for i_cell = 1:number_cells_this_group
    this_cell=cells_this_group(i_cell);
    for i_trial = (trial_counts(i_cell)+1):trial_counts(i_cell+1)
        this_trial_location_IDs=loc_selected(i_cell);
        this_cell=cells_this_group(i_cell);
        this_trial_cell_IDs= this_cell;
        this_trial_power_levels=power_selected(i_cell);
        this_trial_locations=this_neighbourhood.neurons(this_cell).stim_locations.(group_ID).grid(loc_selected(i_cell),:);
            
        
        switch  group_profile.design_func_params.trials_params.stim_design
            case 'Optimal'
                this_trial_power_levels=power_selected(this_cell);
                if rand(1) < mean_gamma(this_cell)
                    this_trial_power_levels=...
                        experiment_setup.prior_info.induced_intensity.fire_stim_threshold./(neurons(i_cell).stim_locations.(group_ID).effect(this_cell,loc_selected(i_cell))...
                        *gain_samples(randsample(1:group_profile.inference_params.MCsamples_for_posterior,1),this_cell));
                    this_trial_power_levels=max(min(power_levels),min(this_trial_power_levels,max(power_levels)));
                end
                
            otherwise
                this_trial_power_levels=randsample(group_profile.design_func_params.trials_params.power_levels,1,true);
        end
            
        experiment_query_this_group.trials(i_trial).location_IDs=this_trial_location_IDs;
        experiment_query_this_group.trials(i_trial).cell_IDs=this_trial_cell_IDs;
        experiment_query_this_group.trials(i_trial).power_levels=this_trial_power_levels;
        experiment_query_this_group.trials(i_trial).locations=this_trial_locations;
        experiment_query_this_group.trials(i_trial).group_ID=group_ID;
        
        
     end
    
end



end
