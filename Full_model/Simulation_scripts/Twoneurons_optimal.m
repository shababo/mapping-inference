
%% Simulate the prior distribution of gain:

% Initialize storage
mean_gamma_current=zeros(n_cell,1);
mean_gain_current=gain_template*ones(n_cell,1);
gamma_path=zeros(n_cell,1);
gain_path=zeros(n_cell,1);
var_gamma_path=zeros(n_cell,1);
var_gain_path=zeros(n_cell,1);

gamma_25_path=zeros(n_cell,1);
gamma_75_path=zeros(n_cell,1);

gain_25_path=zeros(n_cell,1);
gain_75_path=zeros(n_cell,1);


% bounds of the gamma:
gain_bound.up=0.03;gain_bound.low=0.005;
var_alpha_initial=1;var_beta_initial=1.78;
var_alpha_gain_initial=log( (gain_truth+gain_bias- gain_bound.low)./(gain_bound.up-gain_truth));
var_beta_gain_initial=var_gain_prior; % uncertainty of the prior 
% The prioir info serves as the first variational distributions 
variational_params_path.alpha=var_alpha_initial*ones(n_cell,1);
variational_params_path.beta=var_beta_initial*ones(n_cell,1);
variational_params_path.alpha_gain=var_alpha_gain_initial;
variational_params_path.beta_gain=var_beta_gain_initial*ones(n_cell,1);
%% The new selected approach 
n_trials = 0;
n_events=0;
iter=1;

mpp_naive_optimal=cell(0);
    trials_locations_optimal=cell(0);
    trials_powers_optimal=cell(0);
    tic;
fitting_time_optimal=[];

while ((n_trials < trial_max))
    % Select the optimal locations and power levels for each cell:
    
    %---- Select stim locations and power: 
    % draw samples of gains from the posterior distribution 
    
    gain_samples=zeros(n_MC_samples,n_cell);
    for i_cell = 1:n_cell
        v_alpha_gain = variational_params_path.alpha_gain(i_cell,iter);
        v_beta_gain = exp(variational_params_path.beta_gain(i_cell,iter));
        temp=normrnd(v_alpha_gain,v_beta_gain,[n_MC_samples 1]);
        gain_samples(:,i_cell) = exp(temp)./(1+exp(temp))*(gain_bound.up-gain_bound.low) +gain_bound.low; 
    end
    
    firing_prob=zeros(size(pi_target,2),length(power_level),size(pi_target,1));
    for i_loc = 1:size(pi_target,2)
        for k=1:length(power_level)
            for i_cell = 1:n_cell
                 stimulation_received=pi_target(i_cell,i_loc)*power_level(k);
           
                effective_stim= stimulation_received*gain_samples(:,i_cell);
                stim_index=max(1,round(effective_stim*stim_scale));
                prob_collapsed=sum(prob_trace_full(stim_index,:),2);
                firing_prob(i_loc,k,i_cell)=mean(prob_collapsed);
            end
        end
    end
    target_location_optimal=zeros(n_cell,3);
    loc_to_cell_optimal=1:n_cell;
    
    % Select the optimal locations based on firing_prob:
    for i_cell = 1:n_cell
        firing_prob_temp=firing_prob;
        firing_prob_temp(:,:,i_cell)=0;
        firing_prob_difference= firing_prob(:,:,i_cell)-max(firing_prob_temp,[],3);
        [max_value_loc,index_loc] = max(firing_prob_difference);
        % Pick the lowest power if the objectives are not too different from each
        % other
        weighted_max_value_loc = max_value_loc;%./log(log(power_level));
        [~,index_I]=max(weighted_max_value_loc);
        loc_optimal(i_cell)=index_loc(index_I);
        power_optimal(i_cell)=power_level(index_I);
    end

    
    trials_locations=[];
    trials_powers=[];
    
    % draw trials 
    for i_cell = 1:n_cell
        trials_locations( (1:K) +(i_cell-1)*K,1)= loc_optimal(i_cell);
        trials_powers( (1:K) +(i_cell-1)*K,1)= power_optimal(i_cell);
        for k =1:K
            if rand(1) < gamma_path(i_cell,iter)
                trials_powers( (i_cell-1)*K +k,1)=...
                    fire_stim_threshold./(pi_target(i_cell,loc_optimal(i_cell))*gain_samples(k,i_cell));
                trials_powers( (i_cell-1)*K +k,1)=max(min(power_level),min(trials_powers( (i_cell-1)*K +k,1),max(power_level)));
            end
        end
    end
    [~,stim_size] = get_prob_and_size(...
        pi_target,trials_locations,trials_powers,stim_unique,prob_trace);
    
    % Generate mpp given the trials
    [mpp_temp] = draw_samples(...
        trials_locations, trials_powers, pi_target, background_rate,...
        v_th_known, v_reset_known, g_truth, gain_truth,gamma_truth,...
        current_template, funcs, delay_params,stim_threshold,time_max);
    fprintf('There are: %d events;\n',length([mpp_temp.times]))
    mpp_optimal{iter}=mpp_temp;
    trials_locations_optimal{iter}=trials_locations;
    trials_powers_optimal{iter}=trials_powers;
    
    n_trials=n_trials+length(mpp_temp);
    n_events=n_events+length([mpp_temp.times]);
    
    
    %------------------------------------------%
    % Analysis:
    
    % Initialize the path of variational families with infor from previous
    variational_params_path.alpha(:,iter+1)=variational_params_path.alpha(:,iter);
    variational_params_path.beta(:,iter+1)=variational_params_path.beta(:,iter);
    variational_params_path.alpha_gain(:,iter+1)=variational_params_path.alpha_gain(:,iter);
    variational_params_path.beta_gain(:,iter+1)=variational_params_path.beta_gain(:,iter);
    gamma_path(:,iter+1)=gamma_path(:,iter);
    gain_path(:,iter+1)=gain_path(:,iter);
    
    
    variational_params=struct([]);
    for i_cell = 1:n_cell
        variational_params(i_cell).alpha = variational_params_path.alpha(i_cell,iter);
        variational_params(i_cell).beta = variational_params_path.beta(i_cell,iter);
        variational_params(i_cell).alpha_gain = variational_params_path.alpha_gain(i_cell,iter);
        variational_params(i_cell).beta_gain = variational_params_path.beta_gain(i_cell,iter);
    end
    prior_params.pi0= 0.01*ones(length(cell_list),1);
    prior_params.alpha0= [variational_params(:).alpha]';
    prior_params.beta0 = [variational_params(:).beta]';
    prior_params.alpha0_gain= [variational_params(:).alpha_gain]';
    prior_params.beta0_gain =[variational_params(:).beta_gain]';
    
    designs_remained=stim_size;
    mpp_remained=mpp_optimal{iter};
    
    lklh_func=@calculate_likelihood_bernoulli;
    %             lklh_func=@lif_glm_firstevent_loglikelihood_for_VI;
    designs_neighbours=[];        gamma_neighbours=[];         gain_neighbours=[];
    tstart=toc;
    [parameter_history] = fit_working_model_vi_gain(...
        designs_remained, mpp_remained, background_rate, ...
        prob_trace_full,    stim_grid,...
        stim_scale,eff_stim_threshold,gain_bound,...
        variational_params,prior_params,C_threshold,stim_threshold,...
        designs_neighbours,gamma_neighbours,gain_neighbours,...
        S,epsilon,eta_logit,eta_beta,eta_max,maxit,lklh_func);
      tend=toc;
    fitting_time_optimal(iter)=tend-tstart;

    
    [mean_gamma_temp, var_gamma_temp] = calculate_posterior_mean(...
        parameter_history.alpha(:,end),parameter_history.beta(:,end),0,1);
    [mean_gain_temp, var_gain_temp] = calculate_posterior_mean(...
        parameter_history.alpha_gain(:,end),parameter_history.beta_gain(:,end),gain_bound.low,gain_bound.up);
    
    %cell_list(cluster_of_cells{i_cluster})
        variational_params_path.alpha(:,iter+1) = parameter_history.alpha(:,end);
        variational_params_path.beta(:,iter+1) = parameter_history.beta(:,end);
         variational_params_path.alpha_gain(:,iter+1) = parameter_history.alpha_gain(:,end);
        variational_params_path.beta_gain(:,iter+1) = parameter_history.beta_gain(:,end);
   
        var_gamma_path(:,iter+1)=var_gamma_temp;
        gamma_path(:,iter+1)=mean_gamma_temp;
        gain_path(:,iter+1)=mean_gain_temp;
        var_gain_path(:,iter+1)=var_gain_temp;
        
        % for plotting:
%         [gamma_25_path(:,iter+1),gamma_75_path(:,iter+1)] = calculate_posterior_quatiles(...
%             parameter_history.alpha(:,end),parameter_history.beta(:,end),0,1);
%         [gain_25_path(:,iter+1),gain_75_path(:,iter+1)]= calculate_posterior_quatiles(...
%             parameter_history.alpha_gain(:,end),parameter_history.beta_gain(:,end),gain_bound.low,gain_bound.up);
        
       
   iter=iter+1; 
end
%%
params_path_optimal=struct([]);

params_path_optimal(1).alpha=variational_params_path.alpha;
params_path_optimal(1).beta = variational_params_path.beta;
params_path_optimal(1).alpha_gain= variational_params_path.alpha_gain;
params_path_optimal(1).beta_gain= variational_params_path.beta_gain;

var_gamma_path_optimal=var_gamma_path;
var_gain_path_optimal=var_gain_path;
gamma_path_optimal=gamma_path;
gain_path_optimal=gain_path;

n_events_optimal=n_events;