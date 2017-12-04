function [gain_samples_final, gamma_samples_final] = Gibbs_first_spike(mpp, ...
    target_locations, cell_locations,...
    current_template, shape_template,delay_params,linkfunc,...
    stim_threshold, g, background_rate,v_th_known,...
    gain_grid, gain_prior, gamma_grid,gamma_prior,...
    gamma_initial,gain_initial,n_gibbs_sample,...
    n_round_digit)



delay_prob = zeros(delay_params.n_grid+1,1);
if delay_params.type == 1 % normal
    delay_prob = normpdf((0:delay_params.n_grid),delay_params.mean,...
        delay_params.std);
elseif delay_params.type == 2 % gamma
    shape=(delay_params.mean^2)/(delay_params.std^2);
    %scale
    scale = delay_params.mean/shape;
    delay_prob = gampdf((0:delay_params.n_grid),shape,scale);
end
% we approximate the probability with densities
delay_prob = delay_prob/sum(delay_prob);
min_delay = 0;
max_delay = delay_params.n_grid;

n_cell=size(cell_locations,1);
n_trial = length(mpp);
n_grid=length(current_template);


%------------------------------%
% Calculate the size of stimuli
cell_params.locations =  cell_locations;
cell_params.shape_gain = ones(n_cell,1);
cell_template = struct();
cell_template.shape= shape_template;
[pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
    cell_template,target_locations);
stimuli_size=zeros(n_trial,n_cell);
for l = 1:n_trial
    for m = 1:size(mpp(l).locations,2)
        if isnan(mpp(l).locations(m))
        else
            stimuli_size(l,:) = stimuli_size(l,:)+...
                (pi_target(:,mpp(l).locations(m)).*mpp(l).power)';
        end
    end
end

%------------------------------%
% Store the relevant trials for each cell
relevant_trials_per_cell=cell([n_cell 1]);
temp =1:n_trial;
for i_cell = 1:n_cell
    relevant_indicator=stimuli_size(:,i_cell)>stim_threshold;
    relevant_trials_per_cell{i_cell}=temp(relevant_indicator);
end

%------------------------------%
% Calculate the V for unique values of stimuli
stim_unique = unique(stimuli_size);
stim_unique=unique(round(stim_unique,n_round_digit));
num_stim_unique= length(stim_unique);
v_trace=cell([num_stim_unique 1]);
for i_stim = 1:num_stim_unique
    temp_trace = zeros([n_grid 1]);
    stims_temp=current_template*stim_unique(i_stim);
    temp_trace(1)=0;
    for i_t = 2:n_grid
        temp1=reshape(stims_temp(1:(i_t-1)), [i_t-1,1]);
        temp2=reshape(exp( ((1:(i_t-1)) -i_t+1)*g), [1,i_t-1]);
        temp_trace(i_t) = temp2*temp1;
    end
    v_trace{i_stim}=temp_trace;
    fprintf('%d voltage grid done;\n', i_stim);
end
% Evaluate the probability given the gain_grid
prob_trace=cell([num_stim_unique length(gain_grid)]);
n_grid = length(current_template);
v_th_known_one=v_th_known(1); %since they are all the same
for i_stim = 1:num_stim_unique
    v_trace_one=v_trace{i_stim};
    for i_gain = 1:length(gain_grid)
        gain_one=gain_grid(i_gain);
        [prob_first_spike_delayed] = voltage_to_prob(gain_one,  v_trace_one,...
            v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob);
        prob_trace{i_stim,i_gain} =prob_first_spike_delayed;
    end
end

%
prob_by_trials = cell([n_trial 1]);
for i_trial = 1:n_trial
    prob_by_trials{i_trial} =struct('no_spike',[],'spikes',[]);
    prob_by_trials{i_trial}.no_spike = -n_grid*background_rate;
    prob_by_trials{i_trial}.spikes = background_rate*ones([length(mpp(i_trial).times) 1]);
end



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


% Initialize the prob guess

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

%------------------------------%
% Pre-calculate the grid

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