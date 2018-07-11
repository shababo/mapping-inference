function [experiment_query_this_group] = design_undefined(neighbourhood,group_profile,experiment_setup)
% NOTE: we use random design for now, given the adjusted locations 

% Output:
% The experiment_query_this_group struct with fields:
%    experiment_query_this_group.trials(i_trial).location_IDs=this_trial_location_IDs;
%        experiment_query_this_group.trials(i_trial).cell_IDs=this_trial_cell_IDs;
%        experiment_query_this_group.trials(i_trial).power_levels=this_trial_power_levels;
%        experiment_query_this_group.trials(i_trial).locations=this_trial_locations;


%%
disp('designs undef')
group_ID=group_profile.group_ID;
boundary_params=experiment_setup.prior_info.prior_parameters.boundary_params;
radii=group_profile.design_func_params.candidate_grid_params.max_radius;
n_unique_loc=group_profile.design_func_params.trials_params.num_stim_sites;
power_levels=group_profile.design_func_params.trials_params.power_levels;
%i_cells_this_group=find((arrayfun(@(x) strcmp(this_neighbourhood.neurons(:).group_ID,group_ID),this_neighbourhood.neurons(:))));
i_cells_this_group= find(get_group_inds(neighbourhood,group_ID));
number_cells_this_group=length(i_cells_this_group);
number_cells_all= length(neighbourhood.neurons);
loc_counts=zeros(number_cells_this_group,1);

% obtain posterior mean of gamma
batch_ID=neighbourhood.batch_ID;
neurons=neighbourhood.neurons(i_cells_this_group);
properties={'PR'};summary_stat={'mean', 'lower_quantile', 'upper_quantile'};
temp_output=grab_values_from_neurons(batch_ID,neurons,properties,summary_stat);
mean_PR=temp_output.PR.mean;
range_PR=temp_output.PR.upper_quantile-temp_output.PR.lower_quantile;

if group_profile.design_func_params.trials_params.weighted_indicator
    probability_weights = (1-mean_PR)+range_PR;
else
    probability_weights = ones(number_cells_this_group,1);
end


switch group_profile.design_func_params.trials_params.stim_design
    case 'Optimal'
        
%         gain_samples=zeros(group_profile.inference_params.MCsamples_for_posterior,...
%             number_cells_all);
%         for i_cell = 1:number_cells_all
%             alpha_gain=this_neighbourhood.neurons(i_cell).gain_params(end).alpha;
%             beta_gain=this_neighbourhood.neurons(i_cell).gain_params(end).beta;
%             temp=normrnd(alpha_gain,beta_gain,[group_profile.inference_params.MCsamples_for_posterior 1]);
%             gain_samples(:,i_cell) = exp(temp)./(1+exp(temp))*...
%                 range(group_profile.inference_params.bounds.gain)+group_profile.inference_params.bounds.gain(1);
%         end
%         
%         % calculate the firing probability
%         firing_prob=cell(number_cells_this_group,1);
%         
%         for i_cell = 1:number_cells_this_group
%             this_cell=i_cells_this_group(i_cell);
%             candidate_grid=this_neighbourhood.neurons(this_cell).stim_locations.(group_ID);
%             firing_prob{i_cell}=zeros(size(candidate_grid.effect,2),length(group_profile.design_func_params.trials_params.power_levels),...
%                 size(candidate_grid.effect,1));
%             for i_loc = 1:size(candidate_grid.effect,2)
%                 for k=1:length(group_profile.design_func_params.trials_params.power_levels)
%                     for j_cell = 1:size(candidate_grid.effect,1)
%                         if candidate_grid.effect(j_cell,i_loc)>5e-2
%                             stimulation_received=candidate_grid.effect(j_cell,i_loc)*group_profile.design_func_params.trials_params.power_levels(k);
%                             effective_stim= stimulation_received*gain_samples(:,j_cell);
%                             [~, stim_index]=min(abs(effective_stim - experiment_setup.prior_info.induced_intensity.current));
%       
% %                             stim_index=max(1,round(effective_stim*experiment_setup.prior_info.induced_intensity.stim_scale));
%                             prob_collapsed=experiment_setup.prior_info.induced_intensity.prob(stim_index);                            
% %                             sum(experiment_setup.prior_info.induced_intensity.intensity_grid(stim_index,:),2);
%                             firing_prob{i_cell}(i_loc,k,j_cell)=mean(prob_collapsed);
%                         end
%                     end
%                 end
%             end
%         end
%         loc_selected=zeros(number_cells_this_group,1);
%         power_selected=zeros(number_cells_this_group,1);
%         loc_to_cell_selected=1:number_cells_this_group;
%         
%         % Select the optimal locations based on firing_prob:
%         for i_cell = 1:number_cells_this_group
%             this_cell=i_cells_this_group(i_cell);
%             firing_prob_temp=firing_prob{i_cell};
%             firing_prob_temp(:,:,this_cell)=0;
%             firing_prob_difference= firing_prob{i_cell}(:,:,this_cell)-max(firing_prob_temp,[],3);
%             [max_value_loc,index_loc] = max(firing_prob_difference);
%             % Pick the lowest power if the objectives are not too different from each
%             % other
%             weighted_max_value_loc = max_value_loc./log(group_profile.design_func_params.trials_params.power_levels);
%             weighted_max_value_loc( weighted_max_value_loc<0)=0;
%             if max( weighted_max_value_loc)==0
%                 weighted_max_value_loc(:)=1;
%             end
%             index_I = ...
%                 randsample(1:length(weighted_max_value_loc),1,true,weighted_max_value_loc);
%             %         [~,index_I]=max(weighted_max_value_loc);
%             loc_selected(i_cell)=index_loc(index_I);
%             power_selected(i_cell)=group_profile.design_func_params.trials_params.power_levels(index_I);
%         end
    case 'Nuclei'
        loc_selected=ones(number_cells_this_group,1);
        power_selected=zeros(number_cells_this_group,1);
    case 'Random'
        
        loc_selected=[]; 
        loc_mat=zeros(number_cells_this_group,3); 
        loc_list= cell([number_cells_this_group 1]);
        for i_cell = 1:number_cells_this_group
            this_cell=i_cells_this_group(i_cell);
            loc_mat(i_cell,:)= neighbourhood.neurons(this_cell).location +...
             [neighbourhood.neurons(this_cell).posterior_stat(end).shift_x.mean ...
                 neighbourhood.neurons(this_cell).posterior_stat(end).shift_y.mean 0];
             loc_list{i_cell} = zeros(n_unique_loc,3);
             for i_unique_loc = 1:n_unique_loc
                loc_list{i_cell}(i_unique_loc,:)=  loc_mat(i_cell,:)+...
                    [unifrnd(-radii(1), radii(1)) unifrnd(-radii(2), radii(2)) 0];
             end
        end
        % get the adj matrix on loc_mat:
        adj_mat=zeros(number_cells_this_group,number_cells_this_group);
        for i = 1:number_cells_this_group
            for j=1:number_cells_this_group
                rel_loc = loc_mat(i,:)-loc_mat(j,:);
                adj_mat(i,j)=check_in_boundary(rel_loc,boundary_params);
            end
        end
        power_selected=zeros(number_cells_this_group,1);
        
        
    otherwise
        % throw a warning?
end

experiment_query_this_group=struct;
experiment_query_this_group.group_ID=group_ID;
experiment_query_this_group.neighbourhood_ID=neighbourhood.neighbourhood_ID;
experiment_query_this_group.instruction='Compute hologram';
experiment_query_this_group.trials=struct([]);

% DETERMINE NUMBER OF TRIALS
num_trials = min(group_profile.design_func_params.trials_params.trials_per_batch,number_cells_this_group*group_profile.design_func_params.trials_per_cell);


for i_trial = 1:num_trials
    
    this_trial_location_IDs=zeros(1,group_profile.design_func_params.trials_params.spots_per_trial);
    this_trial_cell_IDs=zeros(1,group_profile.design_func_params.trials_params.spots_per_trial);
    this_trial_power_levels=zeros(1,group_profile.design_func_params.trials_params.spots_per_trial);
    this_trial_locations=zeros(group_profile.design_func_params.trials_params.spots_per_trial,3);
    this_trial_adj_power_per_spot = zeros(1,group_profile.design_func_params.trials_params.spots_per_trial);
    
    prob_initial = probability_weights;
    prob_initial = prob_initial./(loc_counts+0.1);
    prob_initial = prob_initial/sum(prob_initial);
    
    
    for i_spot = 1:group_profile.design_func_params.trials_params.spots_per_trial
        try_count = 0;
        loc_found = 0;
        if sum(prob_initial)>0.1
            prob_initial_thresh = prob_initial;
            thresh = median(prob_initial);
            prob_initial_thresh(prob_initial < thresh) = 0;
            while try_count < 10 && ~loc_found
                % location
                switch  group_profile.design_func_params.trials_params.stim_design
                    case 'Random'
                        temp_index = ...
                            randsample(1:number_cells_this_group,1,true,prob_initial_thresh);
                        this_cell = i_cells_this_group(temp_index);
                        temp_loc=randsample(1:n_unique_loc,1);
                       
                    case 'Optimal'
                        temp_index = ...
                            randsample(1:number_cells_this_group,1,true,prob_initial_thresh);
                        this_cell=i_cells_this_group(temp_index);
                        
                        temp_loc =  loc_selected(temp_index);
                    case 'Nuclei'
                        temp_index = ...
                            randsample(1:number_cells_this_group,1,true,prob_initial_thresh);
                        this_cell=i_cells_this_group(temp_index);
                        
                        temp_loc =  loc_selected(temp_index);
                end
                % power
                switch  group_profile.design_func_params.trials_params.stim_design
                    case 'Optimal'
                        this_trial_power_levels(1,i_spot)=power_selected(temp_index);
%                         if rand(1) < mean_gamma(temp_index)
%                             this_trial_power_levels(1,i_spot)=...
%                                experiment_setup.prior_info.induced_intensity.fire_stim_threshold./(neurons(temp_index).stim_locations.(group_ID).effect(this_cell,loc_selected(temp_index))...
%                                 *gain_samples(randsample(1:group_profile.inference_params.MCsamples_for_posterior,1),this_cell));
%                             this_trial_power_levels(1,i_spot)=max(min(power_levels),min(this_trial_power_levels(1,i_spot),max(power_levels)));
%                         end

                    otherwise
                        this_trial_power_levels(1,i_spot)=randsample(group_profile.design_func_params.trials_params.power_levels,1,true);
                        
                end

                if strcmp(experiment_setup.experiment_type,'experiment') || experiment_setup.sim.use_power_calib
                    this_loc = neighbourhood.neurons(this_cell).stim_locations.(group_ID).grid(temp_loc,:);
                    adj_pow_this_loc = round(experiment_setup.exp.ratio_map(round(this_loc(1))+ceil(size(experiment_setup.exp.ratio_map,1)/2),...
                        round(this_loc(2))+ceil(size(experiment_setup.exp.ratio_map,2)/2))*this_trial_power_levels(i_spot));
                    this_trial_total_adj_power_tmp = sum(this_trial_adj_power_per_spot) + adj_pow_this_loc;
                    if ~(this_trial_total_adj_power_tmp > experiment_setup.exp.max_power_ref)
                        loc_found = 1;
                        this_trial_adj_power_per_spot(i_spot) = adj_pow_this_loc;
                    end
                else
                    loc_found = 1;
                    this_trial_adj_power_per_spot(i_spot) = 0;
                end
                try_count = try_count + 1;
            end
        end
        if ~loc_found
            
            this_trial_location_IDs(i_spot)=NaN;
            this_trial_locations(i_spot,:) = NaN;
            this_trial_power_levels(i_spot) = NaN;
            this_trial_cell_IDs(i_spot)=NaN;
            
        else
            
            loc_counts(temp_index)=loc_counts(temp_index)+1;

%             this_trial_location_IDs(i_spot) = temp_loc;
            this_trial_cell_IDs(i_spot) = neighbourhood.neurons(this_cell).cell_ID;   
            this_loc=  loc_list{temp_index}(temp_loc,:);
            this_loc(3)=neighbourhood.center(3);
            this_trial_locations(i_spot,:) = this_loc;
            %this_cell = temp_index;
  
            %prob_initial = subtract_stim_effects(group_ID,temp_index,prob_initial,loc_selected, neurons);
            prob_initial(find(adj_mat(temp_index,:)))=0;
            prob_initial = max(0,prob_initial);

        end
        
    end
    
%     experiment_query_this_group.trials(i_trial).location_IDs = this_trial_location_IDs;
    experiment_query_this_group.trials(i_trial).cell_IDs = this_trial_cell_IDs;
    experiment_query_this_group.trials(i_trial).power_levels = this_trial_power_levels;
    experiment_query_this_group.trials(i_trial).locations = this_trial_locations;
    experiment_query_this_group.trials(i_trial).group_ID = group_ID;
    experiment_query_this_group.trials(i_trial).adj_power_per_spot = this_trial_adj_power_per_spot;
    
end
