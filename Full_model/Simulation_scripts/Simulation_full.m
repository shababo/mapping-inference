%% Online design:
iter=1;
while ((n_trials < trial_max) & (id_continue>0))
    % while not exceeding the set threshold of total trials
    % and there are new cells being excluded
    variational_params=struct([]);
    variational_params(1).alpha = variational_params_path.alpha(:,iter);
    variational_params(1).beta = variational_params_path.beta(:,iter);
    variational_params(1).alpha_gain = variational_params_path.alpha_gain(:,iter);
    variational_params(1).beta_gain = variational_params_path.beta_gain(:,iter);
    gamma_current=gamma_path(:,iter);
    % Conduct random trials
    % On the undefined cells
    mpp_undefined{iter}=[];trials_locations_undefined{iter}=[];trials_powers_undefined{iter}=[];
    if sum(undefined_cells{iter})>0
        cell_list= find(undefined_cells{iter});
        [trials_locations, trials_powers] = random_design(...
            target_locations_selected,power_level,...
            pi_target_selected, inner_normalized_products,single_spot_threshold,...
            variational_params,n_MC_samples,gain_bound,prob_trace_full,gamma_current,  fire_stim_threshold,stim_scale,...
            loc_to_cell,...
            cell_list,n_spots_per_trial,K_undefined,n_replicates);
        [~, stim_size_undefined] = get_prob_and_size(...
            pi_target_selected,trials_locations,trials_powers,stim_unique,prob_trace);
        
        [mpp_temp] = draw_samples(...
            trials_locations, trials_powers, pi_target_selected, background_rate,...
            v_th_known_related, v_reset_known_related, g_related, gain_related,gamma_related,...
            current_template, funcs, delay_params,stim_threshold,time_max);
        % Generate mpp given the trials
        mpp_undefined{iter}=mpp_temp;trials_locations_undefined{iter}=trials_locations;trials_powers_undefined{iter}=trials_powers;
        n_trials=n_trials+size(stim_size_undefined,1);
    end
    
    %-------
    % Conduct trials on group C, the potentially connected cells
    mpp_connected{iter}=[];
    trials_locations_connected{iter}=[];
    trials_powers_connected{iter}=[];
    if sum(potentially_connected_cells{iter})>0
        % Find cells with close to zero gammas
        cell_list= find(potentially_connected_cells{iter});
        [trials_locations,  trials_powers] = random_design(...
                target_locations_nuclei,power_level,...
                pi_target_nuclei, inner_normalized_products,Inf,...
                variational_params,n_MC_samples,gain_bound,prob_trace_full,gamma_current,  fire_stim_threshold,stim_scale,...
                loc_to_cell_nuclei,...
                cell_list,n_spots_per_trial,K_undefined,n_replicates);
        
        [~, stim_size_connected] = get_prob_and_size(...
            pi_target_nuclei,trials_locations,trials_powers,stim_unique,prob_trace);
        % Generate mpp given the trials
        [mpp_temp] = draw_samples(...
            trials_locations, trials_powers, pi_target_nuclei, background_rate,...
            v_th_known_related, v_reset_known_related, g_related, gain_related,gamma_related,...
            current_template,  funcs,    delay_params,stim_threshold,time_max);
        mpp_connected{iter}=mpp_temp;trials_locations_connected{iter}=trials_locations;trials_powers_connected{iter}=trials_powers;
        n_trials=n_trials+size(stim_size_connected,1);
    end
    
    %------------------------------------------%
    % Analysis:
    
    % Initialize the path of variational families with infor from previous
    % iteratin
    variational_params_path.pi(:,iter+1)=var_pi_ini*ones(length(related_cell_list),1);
    variational_params_path.alpha(:,iter+1)=variational_params_path.alpha(:,iter);
    variational_params_path.beta(:,iter+1)=variational_params_path.beta(:,iter);
    variational_params_path.alpha_gain(:,iter+1)=variational_params_path.alpha_gain(:,iter);
    variational_params_path.beta_gain(:,iter+1)=variational_params_path.beta_gain(:,iter);
    
    
    %------------------------------------------------------%
    % Fit VI on Group A: the undefined cells
    mean_gamma_undefined=zeros(length(related_cell_list),1);
    if sum(undefined_cells{iter})>0
        % Find stimulated cells in these trials 
        cell_list= find(sum(stim_size_undefined>10,1)>0);
        % Update variational and prior distribution
        variational_params=struct([]);
        for i_cell_idx = 1:length(cell_list)
            i_cell=cell_list(i_cell_idx);
            variational_params(i_cell_idx).alpha = variational_params_path.alpha(i_cell,iter);
            variational_params(i_cell_idx).beta = variational_params_path.beta(i_cell,iter);
            variational_params(i_cell_idx).alpha_gain = variational_params_path.alpha_gain(i_cell,iter);
            variational_params(i_cell_idx).beta_gain = variational_params_path.beta_gain(i_cell,iter);
        end
        prior_params.pi0= 0.01*ones(length(cell_list),1);
        prior_params.alpha0= [variational_params(:).alpha]';
        prior_params.beta0 = [variational_params(:).beta]';
        prior_params.alpha0_gain= [variational_params(:).alpha_gain]';
        prior_params.beta0_gain =[variational_params(:).beta_gain]';
        
        designs_remained=stim_size_undefined(:,cell_list);
         mpp_remained=mpp_undefined{iter};
           
         designs_neighbours=[];        gamma_neighbours=[];         gain_neighbours=[];
         %     lklh_func=@lif_glm_firstspike_loglikelihood_for_VI;
         lklh_func=@lif_glm_firstevent_loglikelihood_for_VI;
         tstart=toc;
         [parameter_history] = fit_full_model_vi(...
             designs_remained, mpp_remained, background_rate, ...
             prob_trace_full,    stim_grid,...
             stim_scale,eff_stim_threshold,gain_bound,...
             variational_params,prior_params,C_threshold,stim_threshold,...
             designs_neighbours,gamma_neighbours,gain_neighbours,...
             S,epsilon,eta_logit,eta_beta,maxit,lklh_func);
         [mean_gamma_temp, var_gamma_temp] = calculate_posterior_mean(...
             parameter_history.alpha(:,end),parameter_history.beta(:,end),0,1);
         [mean_gain_temp, var_gain_temp] = calculate_posterior_mean(...
             parameter_history.alpha_gain(:,end),parameter_history.beta_gain(:,end),gain_bound.low,gain_bound.up);
         
         %cell_list(cluster_of_cells{i_cluster})
        variational_params_path.alpha(cell_list,iter+1) = parameter_history.alpha(:,end);
        variational_params_path.beta(cell_list,iter+1) = parameter_history.beta(:,end);
        variational_params_path.alpha_gain(cell_list,iter+1) = parameter_history.alpha_gain(:,end);
        variational_params_path.beta_gain(cell_list,iter+1) = parameter_history.beta_gain(:,end);
        
        var_gamma_path(cell_list,iter+1)=var_gamma_temp;
        gamma_path(cell_list,iter+1)=mean_gamma_temp;
        gain_path(cell_list,iter+1)=mean_gain_temp;
        var_gain_path(cell_list,iter+1)=var_gain_temp;
       
        mean_gamma_undefined(cell_list,1)=mean_gamma_temp;
        mean_gamma_current(cell_list)=mean_gamma_temp;
        
    end
    %-------------------------------------------------------------%
    
    %----------------------------------------------%
    % Fit the VI on group C: potentially connected cells
    % This step is different, we shoul fit each neuron seperately if possible
    mean_gamma_connected=zeros(length(related_cell_list),1);
    variance_gamma_connected=ones(length(related_cell_list),1);
    if sum(potentially_connected_cells{iter})>0
        cell_list= find(potentially_connected_cells{iter});
        designs_remained=stim_size_connected(:,cell_list);
        
        % Break the trials into unrelated clusters
        if sum(potentially_connected_cells{iter})>1
            % Use inner product:
            adj_corr= abs( designs_remained'*designs_remained)./size(designs_remained,1);
            adj_corr=1*(adj_corr> (eff_stim_threshold/gain_template/2)^2);
            
            cc_corr=expm(adj_corr);
            cell_cluster_ind=zeros(length(cell_list),1);
            cell_numbers = find(cell_list);
            cluster_id=1;
            for i_cell_idx = 1:length(cell_list)
                
                if cell_cluster_ind(i_cell_idx)==0
                    this_id=cluster_id;
                    cluster_id=cluster_id+1;
                else
                    this_id = cell_cluster_ind(i_cell_idx);
                end
                cell_cluster_ind( find(cc_corr(:,i_cell_idx)))=this_id;
            end
            % Now turn the cell_cluster_ind into list of cells
            n_cluster=max(cell_cluster_ind);
            cluster_of_cells= cell([n_cluster 1]);
            for i_cluster = 1:n_cluster
                cluster_of_cells{i_cluster}=find(cell_cluster_ind==i_cluster);
            end
            
        else
            n_cluster=1;
            cluster_of_cells= cell([n_cluster 1]);
            cluster_of_cells{1}=1;
        end
        
        
        % Now fit the vi model for each of the cluster:
        for i_cluster= 1:n_cluster
%              
            % find the trials that are relevant to thiscell 
            active_trials=find(sum(stim_size_connected(:,cell_list(cluster_of_cells{i_cluster})),2)>stim_threshold);
            neighbour_list= find(sum(stim_size_connected(active_trials,:)>stim_threshold,1)>0);
            
%             neighbour_list=find(sum(cell_neighbours_3d(cell_list(cluster_of_cells{i_cluster}),:),1)>0)';
            variational_params=struct([]);
            for i_cell_idx = 1:length(neighbour_list)
                i_cell=neighbour_list(i_cell_idx);
                 variational_params(i_cell_idx).alpha = variational_params_path.alpha(i_cell,iter);
                variational_params(i_cell_idx).beta = variational_params_path.beta(i_cell,iter);
                variational_params(i_cell_idx).alpha_gain = variational_params_path.alpha_gain(i_cell,iter);
                variational_params(i_cell_idx).beta_gain = variational_params_path.beta_gain(i_cell,iter);
            end
            
            prior_params.pi0= 0.01*ones(length(neighbour_list),1);
        prior_params.alpha0= [variational_params(:).alpha]';
            prior_params.beta0 = [variational_params(:).beta]';
            prior_params.alpha0_gain= [variational_params(:).alpha_gain]';
            prior_params.beta0_gain =[variational_params(:).beta_gain]';
            
            designs_remained=stim_size_connected(active_trials,neighbour_list);
           mpp_remained=mpp_connected{iter}(active_trials);
            
            lklh_func=@lif_glm_firstspike_loglikelihood_for_VI;
%             lklh_func=@lif_glm_firstevent_loglikelihood_for_VI;
            designs_neighbours=[];        gamma_neighbours=[];         gain_neighbours=[];
            [parameter_history] = fit_full_model_vi(...
                designs_remained, mpp_remained, background_rate, ...
                prob_trace_full,    stim_grid,...
                stim_scale,eff_stim_threshold,gain_bound,...
                variational_params,prior_params,C_threshold,stim_threshold,...
                designs_neighbours,gamma_neighbours,gain_neighbours,...
                S,epsilon,eta_logit,eta_beta,maxit,lklh_func);
            
            %cell_list(cluster_of_cells{i_cluster})
        variational_params_path.alpha(neighbour_list,iter+1) = parameter_history.alpha(:,end);
        variational_params_path.beta(neighbour_list,iter+1) = parameter_history.beta(:,end);
         variational_params_path.alpha_gain(neighbour_list,iter+1) = parameter_history.alpha_gain(:,end);
        variational_params_path.beta_gain(neighbour_list,iter+1) = parameter_history.beta_gain(:,end);
       
        [mean_gamma_temp, var_gamma_temp] = calculate_posterior_mean(...
            parameter_history.alpha(:,end),parameter_history.beta(:,end),0,1);
        [mean_gain_temp, ~] = calculate_posterior_mean(...
            parameter_history.alpha_gain(:,end),parameter_history.beta_gain(:,end),gain_bound.low,gain_bound.up);
        
            variance_gamma_connected(neighbour_list)=var_gamma_temp;
            mean_gamma_connected(neighbour_list,1)=mean_gamma_temp;
            mean_gamma_current(neighbour_list)=mean_gamma_temp;
            mean_gain_current(neighbour_list)=mean_gain_temp;
            
            gamma_path(neighbour_list,iter+1)=mean_gamma_temp;
            var_gamma_path(neighbour_list,iter+1)=var_gamma_temp;
            gain_path(neighbour_list,iter+1)=mean_gamma_temp;
            var_gain_path(neighbour_list,iter+1)=var_gamma_temp;
            
            mean_gamma_current(neighbour_list)=mean_gamma_temp;
        
        end
    end
        change_gamma = sqrt(variance_gamma_connected);
        
    %------------------------------------------------------%
    % Moving the cell between groups
    % mean_gamma_undefined & undefined_cells{iter} % A
    % mean_gamma_disconnected & potentially_disconnected_cells{iter} %B
    % mean_gamma_connected & potentially_connected_cells{iter} %C
    
    variational_params=struct([]);
    variational_params(1).alpha = variational_params_path.alpha(:,iter);
    variational_params(1).beta = variational_params_path.beta(:,iter);
    variational_params(1).alpha_gain = variational_params_path.alpha_gain(:,iter);
    variational_params(1).beta_gain = variational_params_path.beta_gain(:,iter);
    
    current_assignments = 1*dead_cells{iter}+2*undefined_cells{iter}+3*potentially_connected_cells{iter}+...
        4*alive_cells{iter};
    moving_thresholds=struct([]);
    moving_thresholds.gamma=[0.1 0.5];
    moving_thresholds.gain_variance=[0.1 0.4];
    
    move_cells(variational_params,current_assignments, gain_bound,moving_thresholds );
    
    %
    iter=iter+1;
    %
    if sum(dead_cells{iter}+alive_cells{iter})==n_cell_this_plane
        id_continue=0;% terminate
    else
        id_continue=1;
    end
    fprintf('Number of trials so far: %d; number of cells confirmed: %d\n',n_trials, sum(dead_cells{iter}+alive_cells{iter}))
end     
  