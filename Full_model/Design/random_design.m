function [ trials_locations,  trials_powers, target_locations_key, pockels_ratio_refs, pockels_ratios] = random_design(...
    target_locations_selected,power_level,...
    pi_target_selected,inner_normalized_products,single_spot_threshold,...
    variational_params,n_MC_samples,gain_bound,prob_trace_full,gamma_current,  fire_stim_threshold,stim_scale,...
    loc_to_cell,...
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
    max_power_ref = varargin{3};
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

pockels_ratio_refs = [];
pockels_ratios = [];

n_remaining_cell=length(remaining_cell_list);
loc_counts=zeros(n_remaining_cell,1);
if n_remaining_cell < single_spot_threshold || n_spots_per_trial==1
    
    n_spots_per_trial=1;
    this_trial_locations=zeros(1,1);
    this_trial_powers=zeros(1,1);
    probability_weights =ones(n_remaining_cell,1);
        
        for i_cell =1:n_remaining_cell % generate trials for each cell
            % Generate random locations on z-plane where the cell sits
            related_locations = find(loc_to_cell == remaining_cell_list(i_cell) );
            trials_temp=[related_locations(1); randsample(related_locations(2:end),K-1,false)];
            
            trials_locations(n_replicates*K*(i_cell-1)+(1:(n_replicates*K)),:)=...
                reshape(repmat(trials_temp,1,n_replicates),[n_replicates*K 1]);
%             randsample(power_level,n_replicates*K,true)'
            trials_powers(n_replicates*K*(i_cell-1)+(1:(n_replicates*K)),:)=...
                randsample(power_level,n_replicates*K,true)';
            
            if use_power_map
                for i = 1:length(trials_temp)
                    for j = 1:n_replicates
                        this_loc = target_locations_selected(trials_temp(i),:);
                        pockels_ratio_refs(end + 1) = round(ratio_map(round(this_loc(1))+ceil(size(ratio_map,1)/2),...
                            round(this_loc(2))+ceil(size(ratio_map,2)/2))*10000)/10000;%/max_power_ref;
                    end
                end
            end
        end
    
else
    % New design for multi-spot stimulation:
    
    % Draw samples from the current posterior distribution of gains 
    gain_samples=zeros(n_MC_samples,length(variational_params.alpha));
    for i_cell = 1:length(variational_params.alpha)
        v_alpha_gain = variational_params.alpha_gain(i_cell);
        v_beta_gain = exp(variational_params.beta_gain(i_cell));
        temp=normrnd(v_alpha_gain,v_beta_gain,[n_MC_samples 1]);
        gain_samples(:,i_cell) = exp(temp)./(1+exp(temp))*(gain_bound.up-gain_bound.low) +gain_bound.low; 
    end
    
    % calculate the firing probability 
    firing_prob=zeros(size(pi_target_selected,2),length(power_level),size(pi_target_selected,1));
    for i_loc = 1:size(pi_target_selected,2)
        for k=1:length(power_level)
            for i_cell = 1:size(pi_target_selected,1)
                if pi_target_selected(i_cell,i_loc)>5e-2
                stimulation_received=pi_target_selected(i_cell,i_loc)*power_level(k);
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
        firing_prob_temp(:,:,i_cell_idx)=0;
        firing_prob_difference= firing_prob(:,:,i_cell)-max(firing_prob_temp,[],3);
        [max_value_loc,index_loc] = max(firing_prob_difference);
        % Pick the lowest power if the objectives are not too different from each
        % other
        weighted_max_value_loc = max_value_loc./log(power_level);
        [~,index_I]=max(weighted_max_value_loc);
        loc_optimal(i_cell_idx)=index_loc(index_I);
        power_optimal(i_cell_idx)=power_level(index_I);
    end
    
    % Design trials using the optimal locations & power 
    n_unique_trials = round(K*n_remaining_cell/n_spots_per_trial);
    trials_locations =zeros(n_unique_trials*n_replicates,n_spots_per_trial);
    trials_powers =zeros(n_unique_trials*n_replicates,n_spots_per_trial);
    n_stim_locations = size(target_locations_selected,1);
    this_trial_locations=zeros(1,n_spots_per_trial);
    
    
    % Set the probability to be inversely proportional to the
    % gamma_estimates, since we want to eliminate more disconnected cells
    probability_weights = ones(n_remaining_cell,1);
    pockels_ratios = zeros(n_unique_trials,n_spots_per_trial);
    for i_trial = 1:n_unique_trials
        this_trial_powers=zeros(1,n_spots_per_trial);
        prob_initial = probability_weights;
        prob_initial = prob_initial./(loc_counts+0.1);
        prob_initial = prob_initial/sum(prob_initial);
        pockels_ratio_refs(end+1) = 0;
        for i_spot = 1:n_spots_per_trial
            
            %            if use_power_map
            %                power_test = pockels_ratio_refs(end) < max_power_ref/power_selected(1)/i_spot; % hack here since all same power...
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
                    temp_index = ...
                        randsample(1:n_remaining_cell,1,true,prob_initial_thresh);
                    temp_loc =  loc_optimal(temp_index);
                    temp_power = power_optimal(temp_index);
                    % with probability proportional to its gamma, jitter the
                    % power:
                    if rand(1) < gamma_current(i_cell)
                        temp_power=...
                        fire_stim_threshold./(pi_target_selected(temp_index,loc_optimal(temp_index))...
                        *gain_samples(randsample(1:n_MC_samples,1),temp_index));
                        temp_power=max(min(power_level),min(temp_power,max(power_level)));
                    end
                    if use_power_map % NOT CODED FOR USE WITH REPLICATING WITHIN THIS FUNCTION!
                        this_loc = target_locations_selected(temp_loc,:);
                        ratio_this_loc = round(ratio_map(round(this_loc(1))+ceil(size(ratio_map,1)/2),...
                            round(this_loc(2))+ceil(size(ratio_map,2)/2))*10000);
                        ratio_this_loc = ratio_this_loc * temp_power;
                        total_adj_power_tmp = sum(pockels_ratios(i_trial,:)/10000) + ratio_this_loc/10000;
                        if ~(total_adj_power_tmp > max_power_ref) 
                            loc_found = 1;
                        end
                    else
                        loc_found = 1;
                        ratio_this_loc=0; % so that the code can run if use_power_map =0;
                        total_adj_power_tmp=0;
                    end
                    try_count = try_count + 1;
                end
            end
            if ~loc_found
                this_trial_locations(1,i_spot)=NaN;
                this_trial_powers(1,i_spot)=0;
            else
                this_trial_locations(1,i_spot)=temp_loc;
                this_trial_powers(1,i_spot)=temp_power;
                loc_counts(temp_index)=loc_counts(temp_index)+1;
                
                pockels_ratios(i_trial,i_spot) = ratio_this_loc;
                pockels_ratio_refs(end) = total_adj_power_tmp/max_power_ref;
                
                
                prob_initial = ...
                    prob_initial - inner_normalized_products(remaining_cell_list,temp_loc)*prob_initial(temp_index);
                prob_initial = max(0,prob_initial);
                
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
            target_locations_key(i,:,k) = target_locations_selected(target_ind,:);
        end
    end
end

end