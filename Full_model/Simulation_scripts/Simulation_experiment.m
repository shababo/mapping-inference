%% Online design:
iter=1;
n_trials=0;
tic
mpp_undefined=struct([]);mpp_connected=struct([]);
computing_time=struct;
computing_time(1).multi=[];
computing_time(1).single=[];
while ((n_trials < trial_max))
    
    variational_params=struct;
    variational_params= variational_params_path(iter,:);
    gamma_current=[parameter_path(iter,:).gamma_mean];
    
    
    % Allocate trials to each groups (by the number of cells in each group)
    % num_trials_per_batch set the maximum number of trials in each batch 
    
    % if total number of trials suggested are large than the threshold
    %   divide the number of trials in proportion
    % otherwise, 
    %   use the suggest number of trials 
    num_suggested_trials = ceil(sum(undefined_cells{iter})*K_undefined/n_spots_per_trial)+...
         sum(connected_cells{iter})*K_connected;
     if num_suggested_trials>num_trials_per_batch
    proportion_undefined=...
        sum(undefined_cells{iter})/(sum(undefined_cells{iter})+multiplier_connected*sum(connected_cells{iter}));
    num_trials_undefined=round(num_trials_per_batch*proportion_undefined);
    num_trials_connected=round(num_trials_per_batch*(1-proportion_undefined));
     else
    num_trials_undefined=ceil(sum(undefined_cells{iter})*K_undefined/n_spots_per_trial);
    num_trials_connected=sum(connected_cells{iter})*K_connected;
         
         
     end
    if num_suggested_trials==0
       break; 
        
    end
    if cell_killing==2
        num_trials_undefined=num_trials_per_batch;
        num_trials_connected=0;
    elseif cell_killing==3
        num_trials_undefined=0;
        num_trials_connected=num_trials_per_batch;
    end
    %---------------------------------------------%
    % Design and conduct trials
    % Undefined cells
    if num_trials_undefined>0
        cell_list=find(undefined_cells{iter});
        
        if cell_killing==2
           cell_list=1:n_cells_this_plane; 
        end
        
        [trials_locations, trials_powers] = random_design(...
            num_trials_undefined, target_locations,power_level,loc_to_cell,related_cell_list,cell_list,...
            pi_target, inner_normalized_products,...
            variational_params,n_MC_samples,gamma_bound,gain_bound,prob_trace_full,...
            fire_stim_threshold,stim_scale, ...
            n_spots_per_trial,[],n_replicates,[],[],[],[],design_type_multi,design_type_single,weighted_design);
        [mpp_temp] = draw_samples(...
            trials_locations, trials_powers, pi_target_truth, background_rate,...
            cell_params(related_cell_list), current_template, funcs, delay_params,stim_threshold,time_max);
        [mpp_temp.batch]=deal(iter);
        if isempty(mpp_undefined)
            mpp_undefined=mpp_temp;
        else
            mpp_undefined(end+(1:+length(mpp_temp)))=mpp_temp;
        end
        n_trials=n_trials+length(mpp_temp);
    end
    
    
    
    %-------
    % Conduct trials on group C, the potentially connected cells
    
    if num_trials_connected>0
        % Find cells with close to zero gammas
        cell_list= find(connected_cells{iter});
          
        if cell_killing==3
           cell_list=1:n_cells_this_plane; 
        end
        [trials_locations,  trials_powers] = random_design(...
                 num_trials_connected, target_locations_nuclei,power_level,loc_to_cell_nuclei,...
                 related_cell_list,cell_list,pi_target_nuclei, [],...
            variational_params,n_MC_samples,gamma_bound,gain_bound,prob_trace_full,...
            fire_stim_threshold,stim_scale, ...
            1,K_connected,n_replicates,[],[],[],[],design_type_multi,design_type_single,weighted_design);
       
       [mpp_temp] = draw_samples(...
            trials_locations, trials_powers, pi_target_nuclei_truth, background_rate,...
            cell_params(related_cell_list), current_template, funcs, delay_params,stim_threshold,time_max);
        [mpp_temp.batch]=deal(iter);
        if isempty(mpp_connected)
            mpp_connected=mpp_temp;
        else
            mpp_connected(end+(1:+length(mpp_temp)))=mpp_temp;
        end
        n_trials=n_trials+length(mpp_temp);
    end
    
    %------------------------------------------%
    % Analysis:
    
    % Initialize the path of variational families with infor from previous
    % iteratin
    variational_params_path(iter+1,:)=variational_params_path(iter,:);
    parameter_path(iter+1,:)=parameter_path(iter,:);
    
tstart=toc;
    %------------------------------------------------------%
    % Fit VI on Group A: the undefined cells
    if num_trials_undefined>0
        % Find stimulated cells in these trials
        
        indicators_remained = find(ismember([mpp_undefined(:).batch],iter-(0:num_trace_back) ));
        mpp_remained=mpp_undefined(indicators_remained);
        trials_locations=reshape([mpp_remained(:).locations],n_spots_per_trial,[])';
        trials_powers=reshape([mpp_remained(:).power],n_spots_per_trial,[])';
        designs_remained = get_stim_size(pi_target,trials_locations,trials_powers);
           
        % include all cells that have been stimulated:
        cell_list= find(sum(designs_remained>stim_threshold,1)>0);
        designs_remained=designs_remained(:,cell_list);
        % Update variational and prior distribution
        variational_params=variational_params_path(iter,cell_list);
        prior_params=variational_params_path(max(iter-num_trace_back,1),cell_list);
        switch model_type
            case 1% working model
                lklh_func=@calculate_loglikelihood_bernoulli;
            case 2% working model with spike & slab
                lklh_func=@calculate_loglikelihood_bernoulli;
            case 3 % full model with first spike
                lklh_func=@lif_glm_firstspike_loglikelihood_for_VI;
            case 4 % full model with first event 
                lklh_func=@lif_glm_firstevent_loglikelihood_for_VI;
            otherwise
        end
        [parameter_history] = fit_VI(...
            designs_remained, mpp_remained, background_rate, ...
            prob_trace_full,stim_scale,minimum_stim_threshold,...
            variational_params,prior_params,gamma_bound,gain_bound,...
            S,epsilon,eta,eta_max,maxit,lklh_func,spike_indicator);
        
        variational_params_path(iter+1,cell_list)=parameter_history(end,:);
        
        [parameter_temp] = calculate_posterior(parameter_history(end,:),gamma_bound,gain_bound,quantiles_prob,spike_indicator);
        parameter_path(iter+1,cell_list)=parameter_temp;
        
        %----------------------------------------%
    end
     tend=toc;
     computing_time(1).multi(iter)=tend-tstart;

    %-------------------------------------------------------------%
    
  
    %----------------------------------------------%
    % Fit the VI on group C: potentially connected cells
    % This step is different, we shoul fit each neuron seperately if possible
     tstart=toc;
    if num_trials_connected>0
        indicators_all = find(ismember([mpp_connected(:).batch],iter-(0:num_trace_back) ));
        mpp_all=mpp_connected(indicators_all);
        trials_locations=reshape([mpp_all(:).locations],1,[])';
        trials_powers=reshape([mpp_all(:).power],1,[])';
        stim_all = get_stim_size(pi_target_nuclei,trials_locations,trials_powers);
        cell_list= find(connected_cells{iter});
        if cell_killing==3
           cell_list=1:n_cells_this_plane; 
        end
      
        [clusters_of_cells] = find_clusters(stim_all, cell_list,stim_threshold);
        % Now fit the vi model for each of the cluster:
        
        for i_cluster= 1:length(clusters_of_cells)
            %
            % find the trials that are relevant to thiscell
            active_trials=find(sum(stim_all(:,cell_list(clusters_of_cells{i_cluster})),2)>stim_threshold);
            neighbour_list= find(sum(stim_all(active_trials,:)>stim_threshold,1)>0);
            variational_params=variational_params_path(iter,neighbour_list);
            prior_params=variational_params_path(max(iter-num_trace_back,1),neighbour_list);
            designs_remained=stim_all(active_trials,neighbour_list);
            mpp_remained=mpp_all(active_trials);
            lklh_func=@lif_glm_firstspike_loglikelihood_for_VI;
            %           lklh_func=@lif_glm_firstevent_loglikelihood_for_VI;
            if model_type == 2
                     spike_indicator=true;  
            else
                spike_indicator=false;
            end
            [parameter_history] = fit_VI(...
                designs_remained, mpp_remained, background_rate, ...
                prob_trace_full,stim_scale,minimum_stim_threshold,...
                variational_params,prior_params,gamma_bound,gain_bound,...
                S,epsilon,eta,eta_max,maxit,lklh_func,spike_indicator);
            

            variational_params_path(iter+1,neighbour_list)=parameter_history(end,:);
            
            [parameter_temp] = calculate_posterior(parameter_history(end,:),gamma_bound,gain_bound,quantiles_prob,spike_indicator);
            parameter_path(iter+1,neighbour_list)=parameter_temp;
            
        end
    end
    tend=toc;
    computing_time(1).single(iter)=tend-tstart;
    
    
    
    % Check the postsynaptic cell type 
    if post_cell_type==1
        [good_cell] = test_celltype(parameter_path,disconnected_params,regular_params,som_params);
        if good_cell
            disconnected_threshold = disconnected_threshold_good;
            connected_threshold =connected_threshold_good;
            quantiles_prob=quantiles_prob_good;
        end
        
    end     
    
    
    [new_list]= move_cells(...
        connected_threshold,disconnected_threshold,single_spot_threshold,multispot_change_threshold,...
        dead_cells{iter},disconnected_cells{iter},undefined_cells{iter},connected_cells{iter},alive_cells{iter},...
        variational_params_path, parameter_path,assignment_type);
    
    disconnected_cells{iter+1}=new_list.disconnected_cells;
    undefined_cells{iter+1}=new_list.undefined_cells;
    connected_cells{iter+1}=new_list.connected_cells;
    alive_cells{iter+1}= new_list.alive_cells;
    dead_cells{iter+1}=new_list.dead_cells;
    % Note: disconnected cells are determined to be dead if they pass the
    % following test:
    
    
    %------------------------------------------------------% 
    % Conduct analysis on group B: potentially disconnected cells
    % - Test for synchrony
    % - Test for small events (OASIS)
    
    if sum(disconnected_cells{iter+1})>0
        %         if small_psc_events == 1
        %
        %         end
        if synchrony_test == 1
            [disconnected_indicators] = test_synchrony(disconnected_cells{iter+1},connected_cells{iter+1}+alive_cells{iter+1},...
                related_cell_list,mpp_connected,mpp_undefined,loc_to_cell,loc_to_cell_nuclei,background_rate,gamma_bound);
            disconnected_to_dead=find(disconnected_indicators & disconnected_cells{iter});
            disconnected_to_connected=find( (~disconnected_indicators) & disconnected_cells{iter});
        else
            disconnected_to_dead=find(disconnected_cells{iter+1});
            disconnected_to_connected=[];
        end
    else
        disconnected_to_dead=[];
        disconnected_to_connected=[];
    end
    %-------------------------------------------------------------%
    dead_cells{iter+1}(disconnected_to_dead)=1;
    disconnected_cells{iter+1}(disconnected_to_dead)=0;
    disconnected_cells{iter+1}(disconnected_to_connected)=0;
    connected_cells{iter+1}(disconnected_to_connected)=1; 
    
    if (sum(dead_cells{iter+1})+sum(alive_cells{iter+1}))  == length(alive_cells{iter+1})
    connected_cells{iter+1}(find(alive_cells{iter+1}))=1;
    end
    
    iter=iter+1;
    
    fprintf('Number of trials: %d; Number of disconnected cells: %d; Number of connected cells: %d\n',n_trials,...
        sum(dead_cells{iter}),sum(alive_cells{iter}))
end
