function [gain_samples_final, gamma_samples_final] = inference_offline(mpp, ...
    target_locations, cell_locations,...
    current_template, shape_template,delay_params,linkfunc,...
    stim_threshold, g, background_rate,v_th_known,...
    gain_grid, gain_prior, gamma_grid,gamma_prior,...
    gamma_initial,gain_initial,n_gibbs_sample,...
    n_round_digit)
% Input: queries, experiment_setup
% Steps:
% 1. Merge all queries: 
% 2. Calculate the stimulation size based o
%     Stim location, cell location (maybe reduce comp cost based on
%     distance?) 
%     Need to leave room for shape estimation
% 3. Initialize gamma and gain, initialize the intensity estimate (for all
% trials)
% 4. Put cells into near-independent clusters
% 5. Run Gibbs sampler
%       - Subtract the effect of one cell from all relevant trials
%       - Draw gain and gamma based on the likelihood 
%       - Add the intensity back 

full_queries
experiment_setup

%% Evaluate the stimulation intensity 
% Need to calculate the stim intensity:



%%

number_cells_this_group=length(i_cell_group_to_nhood);
num_cells_nhood= length(this_neighbourhood.neurons);
prior_info=experiment_setup.prior_info;

%indicators_remained = find(ismember([mpp_undefined(:).batch],iter-(0:num_trace_back) ));
number_of_trials = length(experiment_query_this_group.trials);

% Need a function that graph the mpp from the experiment_query
% note: this_experiment_query contains the group information 
mpp_all=extract_mpp(experiment_query_this_group.trials);
stim_all = get_stim_size(group_ID,experiment_query_this_group.trials,this_neighbourhood);


% include all cells that have been stimulated:
stim_threshold=prior_info.induced_intensity.minimum_stim_threshold/group_profile.inference_params.bounds.gain(2);
%stimulated_cell_list= find(sum(stim_all>stim_threshold,1)>0);

%number_of_stim_cells=length(stimulated_cell_list);
% designs_remained=designs_remained(:,stimulated_cell_list);

%%

n_cell=size(cell_locations,1);
n_trial = length(mpp);
n_grid=length(current_template);

%------------------------------%
% Store the relevant trials for each cell
relevant_trials_per_cell=cell([n_cell 1]);
temp =1:n_trial;
for i_cell = 1:n_cell
    relevant_indicator=stimuli_size(:,i_cell)>stim_threshold;
    relevant_trials_per_cell{i_cell}=temp(relevant_indicator);
end

%% Initialize gamma and gains:


% Initialize the gammas:
if isempty(gamma_initial)
    num_events_trials = zeros(n_trial,1);
    for i_trial = 1:n_trial
        num_events_trials(i_trial) = length(mpp(i_trial).times);
    end
    gamma_initial=zeros(n_cell, 1);
    for i_cell = 1:n_cell
        relevant_trials = relevant_trials_per_cell{i_cell};
        gamma_initial(i_cell)=...
            sum(stimuli_size(relevant_trials,i_cell).*num_events_trials(relevant_trials))...
            /sum(stimuli_size(relevant_trials,i_cell));
        gamma_initial(i_cell) = min(gamma_initial(i_cell),1);
    end
end
% Initialize the gains:
if isempty(gain_initial)
    for i_cell = 1:n_cell
        % the following should be one function:
        relevant_trials =relevant_trials_per_cell{i_cell};
        mpp_relevant=mpp(relevant_trials);
        stim_relevant = stimuli_size(relevant_trials ,i_cell);
        prob_relevant = prob_by_trials(relevant_trials);
        v_th_known_one=v_th_known(i_cell);
        
        
        
        gamma_one=gamma_initial(i_cell);
        if gamma_one == 0
            gain_one_new=randsample(gain_grid,1,true,gain_prior');
        else
            loglklh=zeros([length(gain_grid) 1]);
            for i_grid = 1:length(gain_grid)
                gain_one=gain_grid(i_grid);
                prob_trace_one = prob_trace(:,i_grid);
                [loglklh(i_grid)] = lif_glm_firstspike_loglikelihood_single(gamma_one,...
                    prob_trace_one,stim_unique,stim_relevant,n_round_digit,...
                    prob_relevant,relevant_trials,mpp_relevant,...
                    v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob);
            end
            loggain_post = loglklh;
            loggain_post_shifted=loggain_post-max(loggain_post); % very negative ones become zero
            gain_post=exp(loggain_post_shifted);
            %gain_one_new=randsample(gain_grid,1,true,gain_post);
            [~, gain_one_idx]=max(gain_post);
            gain_one_new = gain_grid(gain_one_idx);
        end
        %         plot(gain_grid,gain_post)
        gain_initial(i_cell)=gain_one_new;
    end
end


% Initialize the prob intensity

for i_cell = 1:n_cell
    relevant_trials =relevant_trials_per_cell{i_cell};
    for i = 1:length(relevant_trials)
        i_trial = relevant_trials(i);
        i_stim = find(stim_unique == round(stimuli_size(i_trial,i_cell),n_round_digit));
        i_gain= find(gain_grid==gain_initial(i_cell));
        prob_first_spike_delayed=prob_trace{i_stim,i_gain};
        
        temp_prob=...
            1-gamma_initial(i_cell)+gamma_initial(i_cell)*(1-sum(prob_first_spike_delayed));
        if temp_prob < 1e-8
            temp_prob=1e-8; % avoid singularity
        end
        prob_by_trials{i_trial}.no_spike =  prob_by_trials{i_trial}.no_spike +...
            log(temp_prob);
        
        if isempty(mpp(i_trial).times)==false
            prob_by_trials{i_trial}.spikes=prob_by_trials{i_trial}.spikes+...
                gamma_initial(i_cell)*prob_first_spike_delayed(round(mpp(i_trial).times));
        end
        
    end
    fprintf('%d prob grid\n', i_cell);
end


%% Put cells into connected components:

[clusters_of_cells] = find_clusters(stim_all, 1:num_cells_nhood, stim_threshold);

%------------------------------%
% Gibbs sampler (a new function)
% iterate through cells


% To further speed up the sampler,
% we can put the cells into separate groups based on
% the stimuli_size matrix
cell_connect_mat = zeros(n_cell,n_cell);
stimuli_thresholded = (stimuli_size > stim_threshold)*1;
cell_connected_mat = stimuli_thresholded'*stimuli_thresholded;
adj_mat=(cell_connected_mat > 0)*1;
% find connected components manually
sparse_adj_mat = sparse(adj_mat);
[n_cc, cc] = graphconncomp(sparse_adj_mat);
%n_cc: number of connected components;
%cc: labels of connected components
cell_list = cell([n_cc 1]);
for i_cc = 1:n_cc
    cell_list{i_cc} = find(cc==i_cc);
end

% Create sparse replicates of the prob_by_trials given the cell_list
prob_by_trials_cc=cell([n_cc 1]);
for i_cc = 1:n_cc
    trial_list = [];
    for i_cell = 1:length(cell_list{i_cc})
        relevant_trials =relevant_trials_per_cell{cell_list{i_cc}(i_cell)};
        trial_list = [trial_list relevant_trials];
    end
    prob_by_trials_cc{i_cc}(unique(trial_list))=prob_by_trials(unique(trial_list));
end

%% Gibbs sampling
gain_samples = cell([n_cc 1]);
gamma_samples = cell([n_cc 1]);

parfor i_cc =1:n_cc
    
    gain_current = gain_initial;
    gamma_current = gamma_initial;
    
    gain_samples{i_cc} = zeros(n_cell,n_gibbs_sample+1);
    gamma_samples{i_cc} = zeros(n_cell,n_gibbs_sample+1);
    gain_samples{i_cc}(:,1)=gain_initial;
    gamma_samples{i_cc}(:,1)=gamma_initial;
    
    for i_sample = 1:n_gibbs_sample
        
        for i_cell_idx = 1:length(cell_list{i_cc})
            % the following should be one function:
            i_cell = cell_list{i_cc}(i_cell_idx);
            relevant_trials =relevant_trials_per_cell{i_cell};
            prob_relevant = prob_by_trials_cc{i_cc}(relevant_trials);
            %v_trace_relevant = v_trace(relevant_trials,i_cell);
            mpp_relevant=mpp(relevant_trials);
            stim_relevant = stimuli_size(relevant_trials ,i_cell);
            
            gain_one_old = gain_current(i_cell);
            gamma_one_old = gamma_current(i_cell);
            v_th_known_one=v_th_known(i_cell);
            
            i_gain= find(gain_grid==gain_one_old);
            prob_trace_one = prob_trace(:,i_gain);
            % Subtract the current cell's effect on the prob_relevant
            [prob_relevant] = subtract_one_cell(gamma_one_old,...
                prob_trace_one,stim_unique,stim_relevant,n_round_digit,...
                prob_relevant,relevant_trials,mpp_relevant,...
                v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob,true);
            % the last one indicate whether to subtract or restore
            
            % Put a function here to calculate the lklh given
            %gain_one and  gamma_one
            gamma_one=gamma_one_old;
            if gamma_one == 0
                gain_one_new=randsample(gain_grid,1,true,gain_prior');
            else
                loglklh=zeros([length(gain_grid) 1]);
                for i_grid = 1:length(gain_grid)
                    gain_one=gain_grid(i_grid);
                    prob_trace_one = prob_trace(:,i_grid);
                    [loglklh(i_grid)] = lif_glm_firstspike_loglikelihood_single(gamma_one,...
                        prob_trace_one,stim_unique,stim_relevant,n_round_digit,...
                        prob_relevant,relevant_trials,mpp_relevant,...
                        v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob);
                end
                loggain_post = loglklh+log(gain_prior');
                loggain_post_shifted=loggain_post-max(loggain_post); % very negative ones become zero
                gain_post=exp(loggain_post_shifted);
                gain_one_new=randsample(gain_grid,1,true,gain_post);
            end
            %         plot(gain_grid,gain_post)
            
            % Calculate the posterior distribution
            % precalculate the prior distribution as
            % draw one from post_gain
            
            i_gain= find(gain_grid==gain_one_new);
            prob_trace_one = prob_trace(:,i_gain);
            loglklh=zeros([length(gamma_grid) 1]);
            for i_grid = 1:length(gamma_grid)
                gamma_one=gamma_grid(i_grid);
                [loglklh(i_grid)] = lif_glm_firstspike_loglikelihood_single(gamma_one,...
                    prob_trace_one,stim_unique,stim_relevant,n_round_digit,...
                    prob_relevant,relevant_trials,mpp_relevant,...
                    v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob);
            end
            loggamma_post = loglklh+log(gamma_prior');
            loggamma_post_shifted=loggamma_post-max(loggamma_post);% very negative ones become zero
            gamma_post=exp(loggamma_post_shifted);
            % draw one from post_gain
            gamma_one_new=randsample(gamma_grid,1,true,gamma_post);
            
            %restore the relevant_prob:
            i_gain= find(gain_grid==gain_one_new);
            prob_trace_one = prob_trace(:,i_gain);
            [prob_relevant] = subtract_one_cell(gamma_one_new,...
                prob_trace_one,stim_unique,stim_relevant,n_round_digit,...
                prob_relevant,relevant_trials,mpp_relevant,...
                v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob,false);
            
            % save it to the main storage
            prob_by_trials_cc{i_cc}(relevant_trials)=prob_relevant;
            gain_current(i_cell)=gain_one_new;
            gamma_current(i_cell)=gamma_one_new;
            %fprintf('Cell %d;', i_cell)
        end
        gain_samples{i_cc}(:,i_sample+1)=gain_current;
        gamma_samples{i_cc}(:,i_sample+1)=gamma_current;
        fprintf('\n The %dth sample drawn in the %d CC;\n', i_sample,i_cc)
    end
    
    
end
gain_samples_final =zeros(n_cell,n_gibbs_sample+1);
gamma_samples_final = zeros(n_cell,n_gibbs_sample+1);
for i_cc = 1:n_cc
    gain_samples_final(cell_list{i_cc},:)=gain_samples{i_cc}(cell_list{i_cc},:);
    gamma_samples_final(cell_list{i_cc},:)=gamma_samples{i_cc}(cell_list{i_cc},:);
end