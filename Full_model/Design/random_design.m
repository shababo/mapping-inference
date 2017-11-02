function [trials_locations,  trials_powers, target_locations_key, pockels_ratio_refs, pockels_ratios] = random_design(...
    num_trials,target_locations,power_level,loc_to_cell,cell_list_map,...
    remaining_cell_list,...
    pi_target,inner_normalized_products,...
    variational_params,n_MC_samples,gamma_bound,gain_bound,prob_trace_full,...
    fire_stim_threshold,stim_scale,...
    n_spots_per_trial,K,n_replicates,...
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

if length(varargin) > 3 && ~isempty(varargin{4})
    do_replicates = varargin{4};
else
    do_replicates = 1;
end

if ~do_replicates
    n_replicates = 1;
    %     K = 1;
end



if length(varargin) > 4 && ~isempty(varargin{5})
    design_type_multi = varargin{5};
else
    design_type_multi = 1;
end


if length(varargin) > 5 && ~isempty(varargin{6})
    design_type_single = varargin{6};
else
    design_type_single = 1;
end


if length(varargin) > 6&& ~isempty(varargin{7})
    weighted_design  = varargin{7};
else
    weighted_design = 1;
end


pockels_ratio_refs = [];
pockels_ratios = [];

n_remaining_cell=length(remaining_cell_list);
loc_counts=zeros(n_remaining_cell,1);
[mean_gamma, ~] = calculate_posterior_mean(...
    [variational_params(remaining_cell_list).alpha],[variational_params(remaining_cell_list).beta],...
    gamma_bound.low,gamma_bound.up);
mean_gamma=mean_gamma.*( 1./(exp([variational_params(remaining_cell_list).p_logit])+1));
if weighted_design 
   probability_weights = 1-mean_gamma;
     
else
    probability_weights =ones(n_remaining_cell,1);
end
    
    
if n_spots_per_trial==1
    trials_locations=[];
    trials_powers=[];
     % Assign the amount of trials per cell based on their weights:
    num_trials_per_cell = round(num_trials*probability_weights/sum(probability_weights));
    num_trials_per_cell(num_trials_per_cell<K)=K;
    if design_type_single==1
        gain_samples=zeros(n_MC_samples,length(variational_params));
        for i_cell = 1:length(variational_params)
            v_alpha_gain = variational_params(i_cell).alpha_gain;
            v_beta_gain = exp(variational_params(i_cell).beta_gain);
            temp=normrnd(v_alpha_gain,v_beta_gain,[n_MC_samples 1]);
            gain_samples(:,i_cell) = exp(temp)./(1+exp(temp))*(gain_bound.up-gain_bound.low) +gain_bound.low;
        end
    end
    for i_cell =1:n_remaining_cell % generate trials for each cell
        % Generate random locations on z-plane where the cell sits
        
        related_locations = find(loc_to_cell == cell_list_map(remaining_cell_list(i_cell)) );
       switch design_type_single
           case 1
               % calculate the firing probability of related locations:
            firing_prob=zeros(length(related_locations),length(power_level),size(pi_target,1));
            for i_loc = 1:length(related_locations)
                for k=1:length(power_level)
                    for idx_cell = 1:size(pi_target,1)
                        if pi_target(idx_cell,related_locations(i_loc))>5e-2
                            stimulation_received=pi_target(idx_cell,related_locations(i_loc))*power_level(k);
                            effective_stim= stimulation_received*gain_samples(:,idx_cell);
                            stim_index=max(1,round(effective_stim*stim_scale));
                            prob_collapsed=sum(prob_trace_full(stim_index,:),2);
                            firing_prob(i_loc,k,idx_cell)=mean(prob_collapsed);
                        end
                    end
                end
            end
            
            % Select the optimal locations based on firing_prob:
            this_cell=remaining_cell_list(i_cell);
            firing_prob_temp=firing_prob;
            firing_prob_temp(:,:,this_cell)=0;
            firing_prob_difference= firing_prob(:,:,this_cell)-max(firing_prob_temp,[],3);
            [max_value_power,index_power] = max(firing_prob_difference');
            max_value_power(max_value_power<0)=0;
            if max(max_value_power(2:end))<1e-3
                max_value_power(:)=1;
            end
            index_loc = [1; ...
                randsample(2:length(max_value_power),num_trials_per_cell(i_cell)-1,true,max_value_power(2:end))'];
            trials_temp=related_locations(index_loc);
            
            trials_locations(end+(1:num_trials_per_cell(i_cell)),:)=reshape(repmat(trials_temp,1,n_replicates),[],1);
            trials_powers(end+(1:num_trials_per_cell(i_cell)),:)=...
                reshape(repmat(power_level(index_power(index_loc)),1,n_replicates),[],1);
      case 0 % random 
               
            trials_temp=[related_locations(1); randsample(related_locations(2:end),...
                num_trials_per_cell(i_cell)-1,true)];
            trials_locations(end+(1:num_trials_per_cell(i_cell)),:)=reshape(repmat(trials_temp,1,n_replicates),[],1);
            trials_powers(end+(1:num_trials_per_cell(i_cell)),:)=randsample(power_level,num_trials_per_cell(i_cell),true)';
           case 2 % nuclei only 
            trials_locations(end+(1:num_trials_per_cell(i_cell)),:)=reshape(repmat(related_locations(1),1,num_trials_per_cell(i_cell)),[],1);
            trials_powers(end+(1:num_trials_per_cell(i_cell)),:)=randsample(power_level,num_trials_per_cell(i_cell),true)';
              
       end
        
        if use_power_map
            for i = 1:length(trials_temp)
                for j = 1:n_replicates
                    this_loc = target_locations(trials_temp(i),:);
                    pockels_ratio_refs(end + 1) = round(ratio_map(round(this_loc(1))+ceil(size(ratio_map,1)/2),...
                        round(this_loc(2))+ceil(size(ratio_map,2)/2))*10000)/10000;
                end
            end
        end
    end
else
    if design_type_multi == 1 % optimal stimulation
        % New design for multi-spot stimulation:
        % Draw samples from the current posterior distribution of gains
        gain_samples=zeros(n_MC_samples,length(variational_params));
        for i_cell = 1:length(variational_params)
            v_alpha_gain = variational_params(i_cell).alpha_gain;
            v_beta_gain = exp(variational_params(i_cell).beta_gain);
            temp=normrnd(v_alpha_gain,v_beta_gain,[n_MC_samples 1]);
            gain_samples(:,i_cell) = exp(temp)./(1+exp(temp))*(gain_bound.up-gain_bound.low) +gain_bound.low;
        end
        % calculate the firing probability
        firing_prob=zeros(size(pi_target,2),length(power_level),size(pi_target,1));
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