% mpp: the marked point process
% event_rates:
% evoked_cell: the list of cells evoked in each trial
% convergence_epsilon: convergence threshold
% mean_background: mean of background events
% sigma_background: standard deviation of background events
% sparsity: indicator for whether to encourage sparsity in the estimates
% gamma_threshold: thresholds for sparse estimates
% sparsity_params: eta (tuning parameter), and threshold (very small value
% to avoid singularity)
function [ ] = Gibbs_first_spike(mpp, ...
    target_locations, cell_locations, current_template, ...
    shape_template,prior_params)
stim_threshold = 10;
gain_initial= 0.03*ones(n_cell,1);
gamma_initial = 0.5*ones(n_cell,1);
g=0.2;
background_rate = 1e-4;
v_th_known=15*ones(n_cell,1);


delay_params.type=2; %1: normal; 2: gamma
delay_params.mean=60;
delay_params.std=30;
delay_params.delayed=true;
delay_params.n_grid=200;

params.invlink = @invlink_sig;
params.dlink = @derlink_sig;
params.link = @link_sig;
params.dinvlink = @derinvlink_sig;
linkfunc = {params.link, params.dlink, params.invlink,params.dinvlink};

 gain_grid=0.001*[1:50];
  gain_prior=normpdf(gain_grid,0.025,0.015);
    gamma_grid= 0.1*[0:10];
  gamma_prior=gamma_grid;
  gamma_prior(1)=0.3;
  gamma_prior(2:end)= (1- gamma_prior(1))/(length(gamma_grid)-1);
  
 
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
% Initialization:
%   for each trial, save the sum log of probability of no spiking
%                        the sum probability of firing at each event
v_trace=cell([n_trial n_cell]);
for i_cell = 1:n_cell
    relevant_trials =relevant_trials_per_cell{i_cell};
    for i = 1:length(relevant_trials)
        i_trial = relevant_trials(i);
        temp_trace = zeros([n_grid 1]);
        stims_temp=current_template*stimuli_size(i_trial,i_cell);
        temp_trace(1)=0;
        for i_t = 2:n_grid
            temp1=reshape(stims_temp(1:(i_t-1)), [i_t-1,1]);
            temp2=reshape(exp( ((1:(i_t-1)) -i_t+1)*g), [1,i_t-1]);
            temp_trace(i_t) = temp2*temp1;
        end
        v_trace{i_trial,i_cell}=temp_trace;
    end
    fprintf('Cell %d voltage grid done;\n', i_cell);
end

% 
prob_by_trials = cell([n_trial 1]);
for i_trial = 1:n_trial
   prob_by_trials{i_trial} =struct('no_spike',[],'spikes',[]);
   prob_by_trials{i_trial}.no_spike = -n_grid*background_rate;
   prob_by_trials{i_trial}.spikes = background_rate*ones([length(mpp(i_trial).times) 1]);
end

for i_cell = 1:n_cell
    relevant_trials =relevant_trials_per_cell{i_cell};
    for i = 1:length(relevant_trials)
        i_trial = relevant_trials(i);        
        prob_firing=zeros([n_grid,1]);
        prob_first_spike=zeros([n_grid,1]);
        prob_first_spike_delayed=zeros([n_grid,1]);
        not_spike_prob=1;
      for i_grid = 2:n_grid
            prob_firing(i_grid)=...
                linkfunc{3}(gain_initial(i_cell)*v_trace{i_trial,i_cell}(i_grid)-...
                v_th_known(i_cell));
            not_spike_prob = not_spike_prob*(1-prob_firing(i_grid -1));
            prob_first_spike(i_grid) =not_spike_prob*prob_firing(i_grid);
        end
        for i_grid = 1:n_grid
            idx_time = max(i_grid-max_delay,1): min(i_grid-min_delay,delay_params.n_grid);
            idx_delay = -( (min(idx_time)-i_grid) : (max(idx_time)-i_grid))+1;
            temp=0;
            for i_time = 1:length(idx_time)
                temp=temp+...
                    prob_first_spike(idx_time(i_time))*delay_prob(idx_delay(i_time));
            end
            prob_first_spike_delayed(i_grid) =temp;
        end
        prob_by_trials{i_trial}.no_spike =  prob_by_trials{i_trial}.no_spike +...
            log(1-gamma_initial(i_cell)+gamma_initial(i_cell)*(1-sum(prob_first_spike_delayed)));
        if isempty(mpp(i_trial).times)==false
              prob_by_trials{i_trial}.spikes=prob_by_trials{i_trial}.spikes+...
                  gamma_initial(i_cell)*prob_first_spike_delayed(round(mpp(i_trial).times));
        end
        
    end
    fprintf('%d prob grid\n', i_cell);
end


%------------------------------%
% Gibbs sampler (a new function) 
% iterate through cells 


gain_current = gain_initial;
gamma_current = gamma_initial;

for i_cell = 1:n_cell
    % the following should be one function:
    
    relevant_trials =relevant_trials_per_cell{i_cell};
    prob_relevant = prob_by_trials(relevant_trials);
    v_trace_relevant = v_trace(relevant_trials,i_cell);
    mpp_relevant=mpp(relevant_trials);
    
    gain_one_old = gain_current(i_cell);
    gamma_one_old = gamma_current(i_cell);
    v_th_known_one=v_th_known(i_cell);
    
    % Subtract the current cell's effect on the prob_relevant
    [prob_relevant] = subtract_one_cell(gain_one_old,gamma_one_old,...
        prob_relevant,v_trace_relevant,relevant_trials,mpp_relevant,...
        v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob,true);
    % the last one indicate whether to subtract or restore
    
    
    % Put a function here to calculate the lklh given
    %gain_one and  gamma_one
    gamma_one=gamma_one_old;
    lklh=zeros([length(gain_grid) 1]);
    for i_grid = 1:length(gain_grid)
        gain_one=gain_grid(i_grid);
        [lklh(i_grid)] = lif_glm_firstspike_loglikelihood_single(gain_one,gamma_one,...
            prob_relevant,v_trace_relevant,relevant_trials,mpp_relevant,...
            v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob);
    end
    
    % Calculate the posterior distribution
    % precalculate the prior distribution as
    gain_post = lklh.*gain_prior';
    % draw one from post_gain
    gain_one_new=randsample(gain_grid,1,true,gain_post);
    
    
    gain_one=gain_one_new;
    lklh=zeros([length(gamma_grid) 1]);
    for i_grid = 1:length(gamma_grid)
        gamma_one=gamma_grid(i_grid);
        [lklh(i_grid)] = lif_glm_firstspike_loglikelihood_single(gain_one,gamma_one,...
            prob_relevant,v_trace_relevant,relevant_trials,mpp_relevant,...
            v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob);
    end
    
    % the draw from the distribution
    gamma_post = lklh.*gamma_prior';
    % draw one from post_gain
    gamma_one_new=randsample(gamma_grid,1,true,gamma_post);
    
    %restore the relevant_prob:
    [prob_relevant] = subtract_one_cell(gain_one_new,gamma_one_new,...
        prob_relevant,v_trace_relevant,relevant_trials,mpp_relevant,...
        v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob,false);
    
    % save it to the main storage
    prob_by_trials(relevant_trials)=prob_relevant;
    gain_current(i_cell)=gain_one_new;
    gamma_current(i_cell)=gamma_one_new;
     fprintf('Draw from cell %d \n', i_cell)
end

end
