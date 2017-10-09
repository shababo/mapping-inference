function data = run_vi_online(data)

params = data.params;
i = data.design.i;

cell_group_list = data.cells_targets.cell_group_list{i};
n_cell_this_plane = length(cell_group_list);

% pi_target_selected = data.cells_targets.pi_target_selected{i};
% inner_normalized_products = data.cells_targets.inner_normalized_products{i};
% target_locations_selected = data.cells_targets.target_locations_selected{i};
% power_selected = data.cells_targets.power_selected{i};
% target_locations_all = data.cells_targets.target_locations_all{i};
cell_neighbours = data.cells_targets.cell_neighbours{i};
% target_locations_nuclei = data.cells_targets.target_locations_nuclei{i};
% power_nuclei = data.cells_targets.power_nuclei{i};
% pi_target_nuclei = data.cells_targets.pi_target_nuclei{i};
% loc_to_cell_nuclei = data.cells_targets.loc_to_cell_nuclei{i};

oasis_data = data.oasis_data;
if ~params.design.do_connected_vi
    full_seq = data.full_seq([data.full_seq.group] ~= 3);
else
    full_seq = data.full_seq;
end
for trial_i = 1:size(oasis_data,1)
    mpp(trial_i).group = full_seq(trial_i).group;
    mpp(trial_i).times = ...
        find(oasis_data(trial_i,...
                params.time.min_time:params.time.max_time),1) + params.time.min_time - 1;
    group_trial_id = full_seq(trial_i).group_target_index;
    switch full_seq(trial_i).group
        case 1
            locations = data.design.trials_locations_undefined{i}{data.design.iter}(group_trial_id,:);
        case 2
            locations = data.design.trials_locations_disconnected{i}{data.design.iter}(group_trial_id,:);
        case 3
            locations = data.design.trials_locations_connected{i}{data.design.iter}(group_trial_id,:);
    end

    mpp(trial_i).locations = locations;
    mpp(trial_i).group_target_index = full_seq(trial_i).group_target_index;
end
assignin('base','mpp',mpp)
data.design.mpp_undefined{i}{data.design.iter} = mpp([mpp.group] == 1);
data.design.mpp_disconnected{i}{data.design.iter} = mpp([mpp.group] == 2);
data.design.mpp_connected{i}{data.design.iter} = mpp([mpp.group] == 3);

%------------------------------------------%
% Transform the data
% no need to record the probabilities all the time..

%data.design.cells_probabilities_undefined;
%         assignin('base','data.design.mpp_undefined{i}',data.design.mpp_undefined{i})
%         assignin('base','data.design.cells_probabilities_undefined',data.design.cells_probabilities_undefined)

% if sum(data.design.undefined_cells{i}{data.design.iter})>0
%     for i_trial = 1:length(data.design.mpp_undefined{i}{data.design.iter})
%         outputs_undefined(i_trial,1)=length(data.design.mpp_undefined{i}{data.design.iter}(i_trial).times);
%     end
%     binary_resp = sort(outputs_undefined > 0);
%     undefined_baseline = mean(binary_resp(1:ceil(length(binary_resp/15))));
%     data.design.n_trials{i}=data.design.n_trials{i}+i_trial;
% end
% if  sum(data.design.potentially_disconnected_cells{i}{data.design.iter})>0
%     %data.design.cells_probabilities_disconnected;
%     for i_trial = 1:length(data.design.mpp_disconnected{i}{data.design.iter})
%         outputs_disconnected(i_trial,1)=length(data.design.mpp_disconnected{i}{data.design.iter}(i_trial).times);
%     end
%     binary_resp = sort(outputs_disconnected > 0);
%     disconnected_baseline = mean(binary_resp(1:ceil(length(binary_resp/15))));
%     data.design.n_trials{i}=data.design.n_trials{i}+i_trial;
% end
% if  sum(data.design.potentially_connected_cells{i}{data.design.iter})>0
%     %data.design.cells_probabilities_disconnected;
% %         for i_trial = 1:size(data.design.stim_size_connected,1)
% %             outputs_connected(i_trial,1)=length(data.design.mpp_connected{i}{data.design.iter}(i_trial).times);
% %         end
%     data.design.n_trials{i}=data.design.n_trials{i}+length(data.design.mpp_connected{i}{data.design.iter});
% end

%------------------------------------------%
% Analysis:

data.design.variational_params_path{i}.pi(:,data.design.iter+1)=params.design.var_pi_ini*ones(n_cell_this_plane,1);
data.design.variational_params_path{i}.alpha(:,data.design.iter+1)=data.design.variational_params_path{i}.alpha(:,data.design.iter);
data.design.variational_params_path{i}.beta(:,data.design.iter+1)=data.design.variational_params_path{i}.beta(:,data.design.iter);
data.design.variational_params_path{i}.alpha_gain(:,data.design.iter+1)=data.design.variational_params_path{i}.alpha_gain(:,data.design.iter);
data.design.variational_params_path{i}.beta_gain(:,data.design.iter+1)=data.design.variational_params_path{i}.beta_gain(:,data.design.iter);


%------------------------------------------------------%
% Fit VI on Group A: the undefined cells
data.design.mean_gamma_undefined=zeros(n_cell_this_plane,1);

if sum(data.design.undefined_cells{i}{data.design.iter})>0
    cell_list= find(data.design.undefined_cells{i}{data.design.iter});
    % Update variational and prior distribution
   variational_params=struct([]);
    for i_cell_idx = 1:length(cell_list)
        i_cell=cell_list(i_cell_idx);
%         variational_params(i_cell_idx).pi = data.design.variational_params_path{i}.pi(i_cell,data.design.iter);
%         variational_params(i_cell_idx).p_logit = log(variational_params(i_cell_idx).pi/(1-variational_params(i_cell_idx).pi));
        variational_params(i_cell_idx).alpha = data.design.variational_params_path{i}.alpha(i_cell,data.design.iter);
        variational_params(i_cell_idx).beta = data.design.variational_params_path{i}.beta(i_cell,data.design.iter);
        variational_params(i_cell_idx).alpha_gain = data.design.variational_params_path{i}.alpha_gain(i_cell,data.design.iter);
        variational_params(i_cell_idx).beta_gain = data.design.variational_params_path{i}.beta_gain(i_cell,data.design.iter);
    end
    prior_params.pi0= 0.01*ones(length(cell_list),1);%[variational_params(:).pi]';
    prior_params.alpha0= [variational_params(:).alpha]';
    prior_params.beta0 = [variational_params(:).beta]';
    prior_params.alpha0_gain= [variational_params(:).alpha_gain]';
    prior_params.beta0_gain =[variational_params(:).beta_gain]';
        
    trial_inds = [data.design.mpp_undefined{i}{data.design.iter}.group_target_index];
    stim_size_undefined = data.design.stim_size_undefined{i}{data.design.iter}(trial_inds,:);
%     cells_probabilities_undefined = cells_probabilities_undefined([data.design.mpp_undefined{i}{data.design.iter}(i_trial).times,:);
    designs_remained = stim_size_undefined(:,cell_list);
    mpp_remained = data.design.mpp_undefined{i}{data.design.iter};
%     active_trials=find(sum(designs_remained,2)>1e-3);
%     designs_remained=designs_remained(active_trials,:);
%     outputs_remained=outputs_undefined(active_trials,:);

    % find neighbours that are not in cell_list:
%     neighbour_list=find(sum(cell_neighbours(cell_list,:),1)>0)';
%     neighbour_list=setdiff(neighbour_list,cell_list);
%     designs_neighbours=cells_probabilities_undefined(active_trials,neighbour_list);
%     gamma_neighbours=data.design.mean_gamma_current{i}(neighbour_list);
    designs_neighbours=[];        gamma_neighbours=[];         gain_neighbours=[];
    lklh_func=@calculate_likelihood_bernoulli;
    % calculate_likelihood_bernoulli for multiple events 
    [parameter_history] = fit_working_model_vi_gain(...
        designs_remained,mpp_remained,params.design.background_rt, ...
        params.prob_trace_full,params.stim_grid,...
        params.stim_scale,params.eff_stim_threshold,params.gain_bound,...
        variational_params,prior_params,params.design.C_threshold,...
        designs_neighbours,gamma_neighbours,gain_neighbors,...
        params.design.S,params.design.epsilon,params.design.eta_logit,...
        params.design.eta_beta,params.design.maxit,lklh_func);
    
    [mean_gamma_temp, var_gamma_temp] = calculate_posterior_mean(...
        parameter_history.alpha(:,end),parameter_history.beta(:,end),0,1);
    [mean_gain_temp, var_gain_temp] = calculate_posterior_mean(...
        parameter_history.alpha_gain(:,end),parameter_history.beta_gain(:,end),gain_bound.low,gain_bound.up);

    % Record the variational parameters
    data.design.variational_params_path.alpha(cell_list,data.design.iter+1) = parameter_history.alpha(:,end);
    data.design.variational_params_path.beta(cell_list,data.design.iter+1) = parameter_history.beta(:,end);
    data.design.variational_params_path.alpha_gain(cell_list,data.design.iter+1) = parameter_history.alpha_gain(:,end);
    data.design.variational_params_path.beta_gain(cell_list,data.design.iter+1) = parameter_history.beta_gain(:,end);

    data.design.var_gamma_path{i}(cell_list,data.design.iter+1)=var_gamma_temp;
    data.design.gamma_path{i}(cell_list,data.design.iter+1)=mean_gamma_temp;
    data.design.gain_path{i}(cell_list,data.design.iter+1)=mean_gain_temp;
    data.design.var_gain_path{i}(cell_list,data.design.iter+1)=var_gain_temp;

    data.design.mean_gamma_undefined(cell_list,1)=mean_gamma_temp;
    data.design.mean_gamma_current{i}(cell_list)=mean_gamma_temp;
    data.design.gamma_path{i}(cell_list,data.design.iter+1)=mean_gamma_temp;

end
%-------------------------------------------------------------%

%----------------------------------------------------------------%
% Fit the VI on Group B: potentially disconnected cells
data.design.mean_gamma_disconnected=ones(n_cell_this_plane,1);
if sum(data.design.potentially_disconnected_cells{i}{data.design.iter})>0
    cell_list= find(data.design.potentially_disconnected_cells{i}{data.design.iter});
    variational_params=struct([]);
    for i_cell_idx = 1:length(cell_list)
        i_cell=cell_list(i_cell_idx);
        variational_params(i_cell_idx).pi = data.design.variational_params_path{i}.pi(i_cell,data.design.iter);
        variational_params(i_cell_idx).p_logit = log(variational_params(i_cell_idx).pi/(1-variational_params(i_cell_idx).pi));
        variational_params(i_cell_idx).alpha = data.design.variational_params_path{i}.alpha(i_cell,data.design.iter);
        variational_params(i_cell_idx).beta = data.design.variational_params_path{i}.beta(i_cell,data.design.iter);
    end

    prior_params.pi0= [variational_params(:).pi]';
    prior_params.alpha0= [variational_params(:).alpha]';
    prior_params.beta0 = [variational_params(:).beta]';
    % Include only the remaining cells
    
    trial_inds = [data.design.mpp_disconnected{i}{data.design.iter}.group_target_index];
    cells_probabilities_disconnected = data.design.cells_probabilities_disconnected{i}{data.design.iter}(trial_inds,:);
    designs_remained=cells_probabilities_disconnected(:,cell_list);
    active_trials=find(sum(designs_remained,2)>1e-3);
    designs_remained=designs_remained(active_trials,:);
    outputs_remained=outputs_disconnected(active_trials,:);

     % find neighbours that are not in cell_list:
    neighbour_list=find(sum(cell_neighbours(cell_list,:),1)>0)';
    neighbour_list=setdiff(neighbour_list,cell_list);
    designs_neighbours=cells_probabilities_disconnected(active_trials,neighbour_list);
    gamma_neighbours=data.design.mean_gamma_current{i}(neighbour_list);

    lklh_func=@calculate_likelihood_bernoulli;
    [parameter_history,~] = fit_working_model_vi(...
        designs_remained,outputs_remained,params.design.background_rt, ...
        variational_params,prior_params,params.design.C_threshold,...
        designs_neighbours,gamma_neighbours,...
        params.design.S,params.design.epsilon,params.design.eta_logit,...
        params.design.eta_beta,params.design.maxit,lklh_func);

    % Record the variational parameters
    data.design.variational_params_path{i}.pi(cell_list,data.design.iter+1) = parameter_history.pi(:,end);
    data.design.variational_params_path{i}.alpha(cell_list,data.design.iter+1) = parameter_history.alpha(:,end);
    data.design.variational_params_path{i}.beta(cell_list,data.design.iter+1) = parameter_history.beta(:,end);

    [mean_gamma_temp, ~] = calculate_posterior_mean(parameter_history.alpha(:,end),parameter_history.beta(:,end),0,1);

    data.design.mean_gamma_disconnected(cell_list,1)=mean_gamma_temp;
    data.design.mean_gamma_current{i}(cell_list)=mean_gamma_temp;
    data.design.gamma_path{i}(cell_list,data.design.iter+1)=mean_gamma_temp;
end
%---------------------------------------------%

%----------------------------------------------%
% Fit the VI on group C: potentially connected cells
% This step is different, we shoul fit each neuron seperately if possible
data.design.mean_gamma_connected=zeros(n_cell_this_plane,1);
data.design.variance_gamma_connected=ones(n_cell_this_plane,1);
lklh_func=@lif_glm_firstspike_loglikelihood_for_VI;
if sum(data.design.potentially_connected_cells{i}{data.design.iter})>0
    cell_list= find(data.design.potentially_connected_cells{i}{data.design.iter});
    
    trial_inds = [data.design.mpp_connected{i}{data.design.iter}.group_target_index];
    stim_size_connected = data.design.stim_size_connected{i}{data.design.iter}(trial_inds,:);
    
    designs_remained=stim_size_connected(:,cell_list);

    % Break the trials into unrelated clusters
    if sum(data.design.potentially_connected_cells{i}{data.design.iter})>1
        % Use inner product:
        adj_corr= abs( designs_remained'*designs_remained)./size(designs_remained,1);
        adj_corr=1*(adj_corr> (params.eff_stim_threshold/params.template_cell.gain_template/2)^2);

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
            cell_cluster_ind(find(cc_corr(:,i_cell_idx)))=this_id;
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

        neighbour_list=find(sum(cell_neighbours(cell_list(cluster_of_cells{i_cluster}),:),1)>0)';

        variational_params=struct([]);
        for i_cell_idx = 1:length(neighbour_list)
            i_cell=neighbour_list(i_cell_idx);
            variational_params(i_cell_idx).pi = data.design.variational_params_path{i}.pi(i_cell,data.design.iter);
            variational_params(i_cell_idx).p_logit = log(variational_params(i_cell_idx).pi/(1-variational_params(i_cell_idx).pi));
            variational_params(i_cell_idx).alpha = data.design.variational_params_path{i}.alpha(i_cell,data.design.iter);
            variational_params(i_cell_idx).beta = data.design.variational_params_path{i}.beta(i_cell,data.design.iter);
            variational_params(i_cell_idx).alpha_gain = data.design.variational_params_path{i}.alpha_gain(i_cell,data.design.iter);
            variational_params(i_cell_idx).beta_gain = data.design.variational_params_path{i}.beta_gain(i_cell,data.design.iter);
        end

        prior_params.pi0= [variational_params(:).pi]';
        prior_params.alpha0= [variational_params(:).alpha]';
        prior_params.beta0 = [variational_params(:).beta]';
        prior_params.alpha0_gain= [variational_params(:).alpha_gain]';
        prior_params.beta0_gain =[variational_params(:).beta_gain]';

        designs_remained=stim_size_connected(:,neighbour_list);
        active_trials=find(sum(designs_remained,2)>params.design.stim_threshold);
        designs_remained=designs_remained(active_trials,:);
        mpp_remained=data.design.mpp_connected{i}{data.design.iter}(active_trials);

%             
%             % find neighbours that are not in cell_list:
%             neighbour_list=find(sum(cell_neighbours(cell_list(cluster_of_cells{i_cluster}),:),1)>0)';
%             neighbour_list=setdiff(neighbour_list,cell_list(cluster_of_cells{i_cluster}));
%             designs_neighbours=data.design.stim_size_connected(active_trials,neighbour_list);
%             gamma_neighbours=data.design.mean_gamma_current{i}(neighbour_list);
%               gain_neighbours=data.design.mean_gain_current{i}(neighbour_list);
%       
        designs_neighbours=[];        gamma_neighbours=[];         gain_neighbours=[];
        [parameter_history] = fit_full_model_vi(...
            designs_remained, mpp_remained, params.bg_rate, ...
            params.template_cell.prob_trace_full,    params.stim_grid,...
            params.stim_scale,params.eff_stim_threshold,params.design.gain_bound,...
            variational_params,prior_params,params.design.C_threshold,params.design.stim_threshold,...
            designs_neighbours,gamma_neighbours,gain_neighbours,...
            params.design.S,params.design.epsilon,params.design.eta_logit,...
            params.design.eta_beta,params.design.maxit,lklh_func);


        %      lklh_func=@calculate_likelihood_bernoulli;
        %     [parameter_history,~] = fit_working_model_vi(...
        %             designs_remained,outputs_remained,background_rt, ...
        %             variational_params,prior_params,C_threshold,...
        %             S,epsilon,eta_logit,eta_beta,maxit,lklh_func);
        %

       %cell_list(cluster_of_cells{i_cluster})
        data.design.variational_params_path{i}.pi(neighbour_list,data.design.iter+1) = parameter_history.pi(:,end);
        data.design.variational_params_path{i}.alpha(neighbour_list,data.design.iter+1) = parameter_history.alpha(:,end);
        data.design.variational_params_path{i}.beta(neighbour_list,data.design.iter+1) = parameter_history.beta(:,end);
        data.design.variational_params_path{i}.alpha_gain(neighbour_list,data.design.iter+1) = parameter_history.alpha_gain(:,end);
        data.design.variational_params_path{i}.beta_gain(neighbour_list,data.design.iter+1) = parameter_history.beta_gain(:,end);


    [mean_gamma_temp, var_gamma_temp] = calculate_posterior_mean(...
        parameter_history.alpha(:,end),parameter_history.beta(:,end),0,1);
    [mean_gain_temp, ~] = calculate_posterior_mean(...
        parameter_history.alpha_gain(:,end),parameter_history.beta_gain(:,end),params.design.gain_bound.low,params.design.gain_bound.up);



        data.design.variance_gamma_connected(neighbour_list)=var_gamma_temp;
        data.design.var_gamma_path(neighbour_list,data.design.iter+1)=var_gamma_temp;
        data.design.mean_gamma_connected(neighbour_list,1)=mean_gamma_temp;
        data.design.mean_gamma_current{i}(neighbour_list)=mean_gamma_temp;
        data.design.mean_gain_current{i}(neighbour_list)=mean_gain_temp;
        data.design.gamma_path{i}(neighbour_list,data.design.iter+1)=mean_gamma_temp;

    end
end
data.design.change_gamma = sqrt(data.design.variance_gamma_connected);


%------------------------------------------------------%

% Debug 
%         etimes1=[mpp_remained(trials_locations==16).times];
%         etimes2=[mpp_remained(trials_locations==14).times];
%         
%         figure(2)
%         scatter(etimes1,0.07*ones(length(etimes1),1),'MarkerFaceColor','r')
%         hold on;
%          line(1:time_max, reshape(probs(:,1,1), [time_max 1]),'Color','r','LineStyle','-')
%          line(1:time_max, reshape(probs(:,1,2), [time_max 1]),'Color','r','LineStyle','--')
%         
%          
%          scatter(etimes2,0.05*ones(length(etimes2),1),'MarkerFaceColor','b')
%          line(1:time_max, reshape(probs(:,2,1), [time_max 1]),'Color','b','LineStyle','-')
%          line(1:time_max, reshape(probs(:,2,2), [time_max 1]),'Color','b','LineStyle','--')
%         
%          
%        hold off;
%        
% %        lines( reshape(probs(:,ii,j)
%        ylim([0 0.1])

%
%        gain_check = mean_gain_temp';
%         gamma_check=mean_gamma_temp;
%        gain_check=gain_truth([14 16])';
%        gamma_check=gamma_truth([14 16]);
%        est_stim_received=designs_remained.*(ones(size(designs_remained,1),1)*gain_check);
%         % estimate the number of events:
% %        est_events = mean_gamma_temp;
% %        for i_cell_temp = 1:length(est_events)
% %            est_events(i_cell_temp)=0;
% %            for i_trial = 1:size(designs_remained,1)
% %                est_events(i_cell_temp)=est_events(i_cell_temp)+stim_to_prob(...
% %                    est_stim_received(i_trial,i_cell_temp),stim_unique,prob_trace);
% %            end
% %        end
%     % 
%     stim_set=unique(est_stim_received,'rows');
%     probs=zeros(time_max,size(stim_set,1),size(stim_set,2));
%     for ii = 1:size(stim_set,1)
%        for j= 1:size(stim_set,2)
%                stim_index=max(1,round(stim_set(ii,j)*stim_scale));
%    probs(:,ii,j)=gamma_check(j)*prob_trace_full(stim_index,:);
%            
%        end
%     end
%             
%             
% 

%

%------------------------------------------------------%
% Moving the cell between groups
% data.design.mean_gamma_undefined & data.design.undefined_cells{i}{data.design.iter} % A
% data.design.mean_gamma_disconnected & data.design.potentially_disconnected_cells{i}{data.design.iter} %B
% data.design.mean_gamma_connected & data.design.potentially_connected_cells{i}{data.design.iter} %C

