addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Load the data set for cell locations
load('./Environments/6_3_s2c3_mpp_and_stim_data.mat')
cell_locations=cell_locs;
n_cell = size(cell_locations,1);
clear mpp;
clear stim_pow;
clear target_inds;
%% Generate the cell parameters

rng(12242,'twister');
background_rate=1e-4;
v_th_known=15*ones([n_cell,1]);
v_reset_known=-1e4*ones([n_cell,1]);
g_truth = 0.02*ones([n_cell,1]);

funcs.invlink = @invlink_sig;%@(resp) log(1 + exp(resp));%@(x) exp(x);
time_max =300;
%gamma_truth=zeros(n_cell_this_plane,1);
gamma_truth = (rand([n_cell 1])<0.1).*(0.7+0.3*rand([n_cell 1]));
% gamma_truth(3)=0.8; gamma_truth(10)=0.8;gamma_truth(1)=0.8;gamma_truth(5)=0.8;
gain_truth=0.015+rand([n_cell 1])*0.01;
%% Preprocessing
%% Estimate the first spike probability given a template cell:
%----------- Delay parameters
load('./Environments/chrome-template-3ms.mat');
downsamp=1;max_time=300;power_level = 30:10:100;
num_power_level=length(power_level);
current_template=template(1:downsamp:max_time);
t_vect= 1:1:max_time;
time_max=max_time;

delay_params.type=2; %1: normal; 2: gamma
delay_params.mean=58; delay_params.std=15;
delay_params.delayed=true; delay_params.n_grid=200;

cell_params.gain_template = 1; % for calculation
cell_params.g=0.02;cell_params.v_th_known=15;
linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};
gain_template=0.02;
stim_scale=4/gain_template;
stim_grid = (1:1000)/stim_scale;
stim_unique=(1:1000)/stim_scale/gain_template;
[prob_trace_full,v_trace_full] = get_first_spike_intensity(...
    linkfunc,...
    current_template,stim_grid,cell_params,delay_params);
prob_trace=sum(prob_trace_full,2);
eff_stim_threshold=stim_grid(min(find(prob_trace>0.01)));
fire_stim_threshold=stim_grid(min(find(prob_trace>0.99)));

%% Put the cells into ten non-overlapping groups by their z-coordinates
n_planes = 10;
z_quantiles = quantile(cell_locations(:,3), (1:(n_planes))*0.1);
cell_group_idx = zeros(n_cell,1);
for i_cell = 1:size(cell_locations,1)
    cell_group_idx(i_cell)= sum(cell_locations(i_cell,3)>z_quantiles)+1;
end
cell_group_list = cell(n_planes,1);
for i_plane = 1:n_planes
    cell_group_list{i_plane} = find(cell_group_idx==i_plane);
end
%% Select one plane as the primary plane
%
this_plane =4;

% Select its neighbours as secondary 
neighbour_plane = [this_plane-1 this_plane+1];
n_cell_this_plane=length(cell_group_list{this_plane});
target_cell_list=struct([]);
target_cell_list(1).primary=cell_group_list{this_plane};
% target_cell_list(1).secondary=[cell_group_list{neighbour_plane(1)}; cell_group_list{neighbour_plane(2)}];
target_cell_list(1).secondary=[];
related_cell_list=[target_cell_list.primary; target_cell_list.secondary]; 


v_th_known_related= v_th_known(related_cell_list);
v_reset_known_related=v_reset_known(related_cell_list);
g_related = g_truth(related_cell_list);

gamma_related = gamma_truth(related_cell_list);
% gamma_truth(3)=0.8; gamma_truth(10)=0.8;gamma_truth(1)=0.8;gamma_truth(5)=0.8;
gain_related=gain_truth(related_cell_list);

%% Precalculate the grid:
load('./Environments/l23_template_cell.mat');
load('./Environments/l23_cells_for_sim.mat');
%l23_average_shape
temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;shape_template=l23_average_shape;

r1=5;r2=10;num_per_grid=12;
grid_type=1;
num_per_grid_dense=16;
stim_threshold=eff_stim_threshold/gain_template;
%
[pi_target_selected, inner_normalized_products,target_locations_selected,...
    target_locations_nuclei, pi_target_nuclei, loc_to_cell_nuclei] = ...
   get_stim_locations(...
    target_cell_list,cell_locations,...
    r1,r2,num_per_grid,num_per_grid_dense,shape_template,...
    grid_type);

loc_to_cell = zeros(size(pi_target_selected,2),1);
for i_cell = 1:length(target_cell_list.primary)
    loc_to_cell( (i_cell-1)*(2*num_per_grid+1)+ (1:(2*num_per_grid+1)))=i_cell;     
end     
%% Designing experiment
%-------------------------------------------------%

% Online design:
%% Parameters in the design stage
% Design parameters
n_spots_per_trial = 4;

n_replicates=1; % number of replicates for each trial 
K_undefined=20; % each cell appears approximately K times
K_disconnected=20; 
K_connected=10; 

single_spot_threshold=9; % switch to single spot stimulation (this can be a function of n_spots_per_trial

trial_max=2000;

% threshold for the group movement:
disconnected_threshold = 0.2;disconnected_confirm_threshold = 0.2;
connected_threshold = 0.5;connected_confirm_threshold = 0.5;
change_threshold=0.01; % for potentially connected cells (since the estimated gamma are important in this case)

% Initialize the five cell groups
undefined_cells= cell(0); undefined_cells{1}=ones(n_cell_this_plane,1);%A
potentially_disconnected_cells= cell(0); potentially_disconnected_cells{1}=zeros(n_cell_this_plane,1);%B
dead_cells= cell(0); dead_cells{1}=zeros(n_cell_this_plane,1);%D
potentially_connected_cells= cell(0); potentially_connected_cells{1}=zeros(n_cell_this_plane,1);%C
alive_cells= cell(0);alive_cells{1}=zeros(n_cell_this_plane,1);%E


% Prior distribution
prior_pi0=0.8;

iter=1;
mpp_undefined=cell(0);trials_locations_undefined=cell(0);trials_powers_undefined=cell(0);

mpp_disconnected=cell(0);trials_locations_disconnected=cell(0);trials_powers_disconnected=cell(0);

mpp_connected=cell(0);trials_locations_connected=cell(0);trials_powers_connected=cell(0);

designs_undefined=[];designs_connected=[];designs_disconnected=[];
outputs_undefined=[];outputs_connected=[];outputs_disconnected=[];


% Initialize the variational family
var_pi_ini=0.01;% not used.
var_alpha_initial=1;var_beta_initial=1.78;
var_alpha_gain_initial=1;var_beta_gain_initial=1.78;

variational_params_path.pi=var_pi_ini*ones(length(related_cell_list),1);
variational_params_path.alpha=var_alpha_initial*ones(length(related_cell_list),1);
variational_params_path.beta=var_beta_initial*ones(length(related_cell_list),1);
variational_params_path.alpha_gain=var_alpha_gain_initial*ones(length(related_cell_list),1);
variational_params_path.beta_gain=var_beta_gain_initial*ones(length(related_cell_list),1);

% Initialize tuning parameters in the VI
C_threshold = 0.01;maxit=1000;
S=200;epsilon=0.01;eta_logit=0;
background_rt=background_rate*time_max; % raw probability of firing within a trial

eta_beta=0.05;

% Whether to output plots during the experiment
visualized = 0;

n_trials=0;


gamma_estimates = 0.5*ones(n_cell_this_plane,1);% for drawing samples (not really used)
prob_weight=0;

id_continue=1;% an indicator


% lklh_func=@calculate_likelihood_sum_bernoulli; % likelihood function is
% specificed when fitting the working model

stim_threshold = 10;

% bounds of the gamma:
gain_bound.up=0.03;gain_bound.low=0.01;

% Used in the random designs
id_notconnected=false;
loc_to_cell = 1:size(target_locations_selected,1);
connected=true;
%  loc_to_cell_nuclei is from get_stim_locations

% Initialize storage
mean_gamma_current=zeros(length(related_cell_list),1);
mean_gain_current=gain_template*ones(length(related_cell_list),1);
gamma_path=zeros(length(related_cell_list),1);
var_gamma_path=zeros(length(related_cell_list),1);
gain_path=zeros(length(related_cell_list),1);
var_gain_path=zeros(length(related_cell_list),1);

change_threshold = 0.05;
prob_lower_bound = 0.01;
stim_threshold=10;

%% Initialize the posterior/variational distributions
% bounds of the gamma:
gain_bound.up=0.03;gain_bound.low=0.005;
var_alpha_initial=1;var_beta_initial=1.78;

var_alpha_gain_initial=log( (0.02 - gain_bound.low)./(gain_bound.up-0.02))*ones(length(related_cell_list),1);
var_beta_gain_initial=0.5; % uncertainty of the prior 

% The prioir info serves as the first variational distributions 
variational_params_path.alpha=var_alpha_initial*ones(length(related_cell_list),1);
variational_params_path.beta=var_beta_initial*ones(length(related_cell_list),1);
variational_params_path.alpha_gain=var_alpha_gain_initial;
variational_params_path.beta_gain=var_beta_gain_initial*ones(length(related_cell_list),1);

n_MC_samples=50;
%% Online design:
iter=1;
% while ((n_trials < trial_max) & (id_continue>0))
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
    % Conduct trials on group B, the potentially disconnected cells
%     mpp_disconnected{iter}=[];trials_locations_disconnected{iter}=[];trials_powers_disconnected{iter}=[];
%     if sum(potentially_disconnected_cells{iter})>0
%         % Find cells with close to zero gammas
%         cell_list= find(potentially_disconnected_cells{iter});
%         if length(cell_list) > single_spot_threshold
%             [trials_locations, trials_powers] = random_design(...
%                 target_locations_selected,power_level,...
%                 pi_target_selected, inner_normalized_products,single_spot_threshold,...
%                 variational_params,n_MC_samples,gain_bound,prob_trace_full,gamma_current,  fire_stim_threshold,stim_scale,...
%                 loc_to_cell,...
%                 cell_list,n_spots_per_trial,K_undefined,n_replicates);
%             [~, stim_size_disconnected] = get_prob_and_size(...
%                 pi_target_selected,trials_locations,trials_powers,stim_unique,prob_trace);
%              [mpp_temp] = draw_samples(...
%             trials_locations, trials_powers, pi_target_selected, background_rate,...
%             v_th_known_related, v_reset_known_related, g_related, gain_related,gamma_related,...
%             current_template, funcs, delay_params,stim_threshold,time_max);
%         else
%             [trials_locations,  trials_powers] = random_design(...
%                 target_locations_nuclei,power_level,...
%                 pi_target_nuclei, inner_normalized_products,single_spot_threshold,...
%                 variational_params,n_MC_samples,gain_bound,prob_trace_full,gamma_current,  fire_stim_threshold,stim_scale,...
%                 loc_to_cell_nuclei,...
%                 cell_list,n_spots_per_trial,K_undefined,n_replicates);
%             [~, stim_size_disconnected] = get_prob_and_size(...
%                 pi_target_nuclei,trials_locations,trials_powers,stim_unique,prob_trace);
%             [mpp_temp] = draw_samples(...
%                 trials_locations, trials_powers, pi_target_nuclei, background_rate,...
%                 v_th_known_related, v_reset_known_related, g_related, gain_related,gamma_related,...
%                 current_template,  funcs,    delay_params,stim_threshold,time_max);
%         end
%         mpp_disconnected{iter}=mpp_temp;trials_locations_disconnected{iter}=trials_locations;trials_powers_disconnected{iter}=trials_powers;
%         n_trials=n_trials+size(stim_size_disconnected,1);
%     end
    
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
        %[cells_probabilities_connected, ~] = get_prob_and_size(...
        %    pi_target_nuclei,trials_locations,trials_powers,...
        %    stim_unique,prob_trace);
        [~, stim_size_connected] = get_prob_and_size(...
            pi_target_nuclei,trials_locations,trials_powers,stim_unique,prob_trace);
        % Generate mpp given the trials
        [mpp_temp] = draw_samples(...
            trials_locations, trials_powers, pi_target_nuclei, background_rate,...
            v_th_known_related, v_reset_known_related, g_related, gain_related,gamma_related,...
            current_template,  funcs,    delay_params,stim_threshold,time_max);
        mpp_connected{iter}=mpp_temp;
        trials_locations_connected{iter}=trials_locations;
        trials_powers_connected{iter}=trials_powers;
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
%         cell_list= find(undefined_cells{iter});
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
           
        lklh_func=@calculate_likelihood_bernoulli;
        designs_neighbours=[];        gamma_neighbours=[];         gain_neighbours=[];
        [parameter_history] = fit_working_model_vi_gain(...
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
    
    %----------------------------------------------------------------%
    % Fit the VI on Group B: potentially disconnected cells
%     mean_gamma_disconnected=ones(length(related_cell_list),1);
%     
%     if sum(potentially_disconnected_cells{iter})>0
%         % Find stimulated cells in these trials 
%        cell_list= find(sum(stim_size_disconnected>10,1)>0); 
%        % Update variational and prior distribution
%         variational_params=struct([]);
%         for i_cell_idx = 1:length(cell_list)
%             i_cell=cell_list(i_cell_idx);
%             variational_params(i_cell_idx).alpha = variational_params_path.alpha(i_cell,iter);
%             variational_params(i_cell_idx).beta = variational_params_path.beta(i_cell,iter);
%             variational_params(i_cell_idx).alpha_gain = variational_params_path.alpha_gain(i_cell,iter);
%             variational_params(i_cell_idx).beta_gain = variational_params_path.beta_gain(i_cell,iter);
%         end
%         prior_params.pi0= 0.01*ones(length(cell_list),1);
%         prior_params.alpha0= [variational_params(:).alpha]';
%         prior_params.beta0 = [variational_params(:).beta]';
%         prior_params.alpha0_gain= [variational_params(:).alpha_gain]';
%         prior_params.beta0_gain =[variational_params(:).beta_gain]';
%         
%         designs_remained=stim_size_disconnected(:,cell_list);
%          mpp_remained=mpp_disconnected{iter};
%            
%         lklh_func=@calculate_likelihood_bernoulli;
%         designs_neighbours=[];        gamma_neighbours=[];         gain_neighbours=[];
%         [parameter_history] = fit_working_model_vi_gain(...
%             designs_remained, mpp_remained, background_rate, ...
%             prob_trace_full,    stim_grid,...
%             stim_scale,eff_stim_threshold,gain_bound,...
%             variational_params,prior_params,C_threshold,stim_threshold,...
%             designs_neighbours,gamma_neighbours,gain_neighbours,...
%             S,epsilon,eta_logit,eta_beta,maxit,lklh_func);
%        [mean_gamma_temp, var_gamma_temp] = calculate_posterior_mean(...
%             parameter_history.alpha(:,end),parameter_history.beta(:,end),0,1);
%         [mean_gain_temp, var_gain_temp] = calculate_posterior_mean(...
%             parameter_history.alpha_gain(:,end),parameter_history.beta_gain(:,end),gain_bound.low,gain_bound.up);
%         
%         %cell_list(cluster_of_cells{i_cluster})
%         variational_params_path.alpha(cell_list,iter+1) = parameter_history.alpha(:,end);
%         variational_params_path.beta(cell_list,iter+1) = parameter_history.beta(:,end);
%         variational_params_path.alpha_gain(cell_list,iter+1) = parameter_history.alpha_gain(:,end);
%         variational_params_path.beta_gain(cell_list,iter+1) = parameter_history.beta_gain(:,end);
%         
%       
%        var_gamma_path(cell_list,iter+1)=var_gamma_temp;
%         gamma_path(cell_list,iter+1)=mean_gamma_temp;
%         gain_path(cell_list,iter+1)=mean_gain_temp;
%         var_gain_path(cell_list,iter+1)=var_gain_temp;
%        
%         mean_gamma_disconnected(cell_list,1)=mean_gamma_temp;
%         mean_gamma_current(cell_list)=mean_gamma_temp;
%     end
    %---------------------------------------------%
    
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
%           lklh_func=@lif_glm_firstevent_loglikelihood_for_VI;
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
        [mean_gain_temp, var_gain_temp] = calculate_posterior_mean(...
            parameter_history.alpha_gain(:,end),parameter_history.beta_gain(:,end),gain_bound.low,gain_bound.up);
        
            variance_gamma_connected(neighbour_list)=var_gamma_temp;
            mean_gamma_connected(neighbour_list,1)=mean_gamma_temp;
            mean_gamma_current(neighbour_list)=mean_gamma_temp;
            mean_gain_current(neighbour_list)=mean_gain_temp;
            
            gamma_path(neighbour_list,iter+1)=mean_gamma_temp;
            var_gamma_path(neighbour_list,iter+1)=var_gamma_temp;
            gain_path(neighbour_list,iter+1)=mean_gamma_temp;
            var_gain_path(neighbour_list,iter+1)=var_gain_temp;
        end
    end
%     if vf_type == 1
%     change_gamma =abs(gamma_path(:,iter+1)-gamma_path(:,iter));
%     elseif vf_type == 2
        change_gamma = sqrt(variance_gamma_connected);
%     end


    %------------------------------------------------------%
    % Moving the cell between groups
    % mean_gamma_undefined & undefined_cells{iter} % A
    % mean_gamma_disconnected & potentially_disconnected_cells{iter} %B
    % mean_gamma_connected & potentially_connected_cells{iter} %C
    
    undefined_to_disconnected = intersect(find(mean_gamma_undefined<disconnected_threshold),find( undefined_cells{iter}));
    undefined_to_connected = intersect(find(mean_gamma_undefined>connected_threshold),find( undefined_cells{iter}));
    
    % cells move together with their neighbours
%     undefined_to_disconnected=find(sum(cell_neighbours(undefined_to_disconnected,:),1)>0)';
%     undefined_to_connected =find(sum(cell_neighbours(undefined_to_connected,:),1)>0);
    % if there are conflicts, move them to the potentially connected cells
%     undefined_to_disconnected=setdiff(undefined_to_disconnected,undefined_to_connected);
    
    disconnected_to_undefined = intersect(find(mean_gamma_disconnected>disconnected_confirm_threshold),...
        find(potentially_disconnected_cells{iter}));
    disconnected_to_dead = intersect(find(mean_gamma_disconnected<disconnected_confirm_threshold),...
        find(potentially_disconnected_cells{iter}));
    
%     disconnected_to_undefined=find(sum(cell_neighbours(disconnected_to_undefined,:),1)>0);
    % if there are conflicts, move them to the potentially connected cells
%     disconnected_to_dead=setdiff(disconnected_to_dead,disconnected_to_undefined);
    
    
    connected_to_dead = intersect(find(mean_gamma_connected<disconnected_confirm_threshold),...
        find(potentially_connected_cells{iter}));
    connected_to_alive = intersect(find(mean_gamma_connected>connected_confirm_threshold),...
        find(potentially_connected_cells{iter}));
    connected_to_alive = intersect(find(change_gamma<change_threshold),...
        connected_to_alive);
    
    % Update the cell lists:
    undefined_cells{iter+1}=undefined_cells{iter};
    undefined_cells{iter+1}(undefined_to_disconnected)=0;undefined_cells{iter+1}(undefined_to_connected)=0;
    undefined_cells{iter+1}(disconnected_to_undefined)=1;
    
    potentially_disconnected_cells{iter+1}=potentially_disconnected_cells{iter};
    potentially_disconnected_cells{iter+1}(disconnected_to_dead)=0;potentially_disconnected_cells{iter+1}(disconnected_to_undefined)=0;
    potentially_disconnected_cells{iter+1}(undefined_to_disconnected)=1;
    
    potentially_connected_cells{iter+1}=potentially_connected_cells{iter};
    potentially_connected_cells{iter+1}(connected_to_dead)=0;potentially_connected_cells{iter+1}(connected_to_alive)=0;
    potentially_connected_cells{iter+1}(undefined_to_connected)=1;
    
    dead_cells{iter+1}=dead_cells{iter};
    dead_cells{iter+1}(disconnected_to_dead)=1;dead_cells{iter+1}(connected_to_dead)=1;
    
    alive_cells{iter+1}=alive_cells{iter};
    alive_cells{iter+1}(connected_to_alive)=1;
    
    %
    iter=iter+1;
    %
    if sum(dead_cells{iter}+alive_cells{iter})==n_cell_this_plane
        id_continue=0;% terminate
    else
        id_continue=1;
    end
    fprintf('Number of trials so far: %d; number of cells confirmed: %d\n',n_trials, sum(dead_cells{iter}+alive_cells{iter}))
% end
    %%
     find(gamma_truth(target_cell_list.primary)>0)
     find(alive_cells{iter})
     find(potentially_connected_cells{iter})
     find(undefined_cells{iter})
     
  