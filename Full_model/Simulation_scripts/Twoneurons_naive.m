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

%% All random approach
n_trials = 0;
n_events=0;
iter=1;
mpp_naive=cell(0);
    trials_locations_naive=cell(0);
    trials_powers_naive=cell(0);
    tic;
fitting_time_naive=[];

while ((n_trials < trial_max))
    
    gamma_estimates = 0.5*ones(length(cell_list),1);% for drawing samples...
    [trials_locations, trials_powers] = random_design_v2(target_locations,...
        power_selected,power_sd,...
        inner_normalized_products,1,gamma_estimates,0,...
        true, loc_to_cell,cell_list,n_spots_per_trial,K,n_replicates,...
        power_level);
    [~,stim_size] = get_prob_and_size(...
        pi_target,trials_locations,trials_powers,stim_unique,prob_trace);
    
    % Generate mpp given the trials
    [mpp_temp] = draw_samples(...
        trials_locations, trials_powers, pi_target, background_rate,...
        v_th_known, v_reset_known, g_truth, gain_truth,gamma_truth,...
        current_template, funcs, delay_params,stim_threshold,time_max);
    mpp_naive{iter}=mpp_temp;
    trials_locations_naive{iter}=trials_locations;
    trials_powers_naive{iter}=trials_powers;
     fprintf('There are: %d events;\n',length([mpp_temp.times]))
    n_events=n_events+length([mpp_temp.times]);
    n_trials=n_trials+length(mpp_temp);
    
    
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
    mpp_remained=mpp_naive{iter};
    
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
        S,epsilon,eta_logit,eta_beta,maxit,lklh_func);
    tend=toc;
    fitting_time_naive(iter)=tend-tstart;

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
        
         %for plotting:
%         [gamma_25_path(:,iter+1),gamma_75_path(:,iter+1)] = calculate_posterior_quatiles(...
%             parameter_history.alpha(:,end),parameter_history.beta(:,end),0,1);
%         [gain_25_path(:,iter+1),gain_75_path(:,iter+1)]= calculate_posterior_quatiles(...
%             parameter_history.alpha_gain(:,end),parameter_history.beta_gain(:,end),gain_bound.low,gain_bound.up);
        
   iter=iter+1; 
end

%%
params_path_naive=struct([]);

params_path_naive(1).alpha=variational_params_path.alpha;
params_path_naive(1).beta = variational_params_path.beta;
params_path_naive(1).alpha_gain= variational_params_path.alpha_gain;
params_path_naive(1).beta_gain= variational_params_path.beta_gain;

var_gamma_path_naive=var_gamma_path;
var_gain_path_naive=var_gain_path;
gamma_path_naive=gamma_path;
gain_path_naive=gain_path;

n_events_naive=n_events;