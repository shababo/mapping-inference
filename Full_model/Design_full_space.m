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
%--------------------------------------------------%
%% Estimate the first spike probability given a template cell:
%----------- Delay parameters
delay_params.type=2; %1: normal; 2: gamma
delay_params.mean=58; delay_params.std=15;
delay_params.delayed=true; delay_params.n_grid=200;

%----------- Load the current template
load('./Environments/chrome-template-3ms.mat');
downsamp=1;max_time=300;power_level = [50 75 100];
num_power_level=length(power_level);
current_template=template(1:downsamp:max_time);
t_vect= 1:1:max_time;
gain_template=0.02;
cell_params.g=0.02;cell_params.v_th_known=15;cell_params.gain_template = gain_template;
stim_unique = (1:1000)/10;
linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};
[prob_trace]=get_firing_probability(...
    linkfunc,current_template,stim_unique,cell_params,delay_params);
%% Another precalculation
%------------------------------%
% Calculate the firing intensity
stim_scale=200;
stim_grid = (1:1000)/stim_scale;
cell_params.g=0.02;
[prob_trace_full,v_trace_full] = get_first_spike_intensity(...
    linkfunc,...
    current_template,stim_grid,cell_params,delay_params);
eff_stim_threshold=stim_grid(min(find(sum(prob_trace_full,2)>1e-1)));
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
target_cell_list(1).secondary=[cell_group_list{neighbour_plane(1)}; cell_group_list{neighbour_plane(1)}];
related_cell_list=[target_cell_list.primary; target_cell_list.secondary]; 


v_th_known_related= v_th_known(related_cell_list);
v_reset_known_related=v_reset_known(related_cell_list);
g_related = g_truth(related_cell_list);

gamma_related = gamma_truth(related_cell_list);
% gamma_truth(3)=0.8; gamma_truth(10)=0.8;gamma_truth(1)=0.8;gamma_truth(5)=0.8;
gain_related=gain_truth(related_cell_list);

%% Select the stimulation locations
% Load the shape template
load('./Environments/l23_template_cell.mat');
load('./Environments/l23_cells_for_sim.mat');
%l23_average_shape
temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;shape_template=l23_average_shape;
r1=5;r2=10;r3=15;num_per_grid=12;
num_per_grid_dense=16;
stim_threshold=eff_stim_threshold/gain_template;
grid_type = 2; % 1: circles; 2: lines;
[pi_target_selected, inner_normalized_products,target_locations_selected,power_selected,...
    target_locations_all,cell_neighbours,...
    target_locations_nuclei, power_nuclei,pi_target_nuclei, loc_to_cell_nuclei] = ...
    get_stim_locations(...
    target_cell_list,cell_locations,power_level,...
    r1,r2,r3,num_per_grid,num_per_grid_dense,shape_template,...
    stim_unique,prob_trace,stim_threshold,...
    grid_type);
%num_z_grid=8;% 6 values on z-plane (1 N N N N 1)
%%

%% End of preprocessing
%-------------------------------------------------%

%% Designing experiment
%-------------------------------------------------%

% Online design:
%% Parameters in the design stage

% Design parameters
n_spots_per_trial = 4;

% Need to run sims to check how these parameters affect the results

n_replicates=1; % conduct two replicates for each trial
K_undefined=5; % each cell appears approximately 10*2 times
K_disconnected=5; % each cell appears approximately 10*2 times
K_connected=10; % each cell appears approximately 10*2 times


single_spot_threshold=15; % switch to single spot stimulation (this can be a function of n_spots_per_trial

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

change_threshold = 0.05;
prob_lower_bound = 0.01;
stim_threshold=10;
%% Online design:
while ((n_trials < trial_max) & (id_continue>0))
    % while not exceeding the set threshold of total trials
    % and there are new cells being excluded
    
    % Conduct random trials
    % On the undefined cells
    mpp_undefined{iter}=[];trials_locations_undefined{iter}=[];trials_powers_undefined{iter}=[];
    if sum(undefined_cells{iter})>0
        cell_list= find(undefined_cells{iter});
        gamma_estimates = 0.5*ones(length(cell_list),1);% for drawing samples...
        [trials_locations, trials_powers] = random_design(target_locations_selected,power_selected,...
            inner_normalized_products,single_spot_threshold,gamma_estimates,prob_weight,...
            id_notconnected, loc_to_cell,cell_list,n_spots_per_trial,K_undefined,n_replicates);
        [cells_probabilities_undefined, ~] = get_prob_and_size(...
            pi_target_selected,trials_locations,trials_powers,stim_unique,prob_trace);
        
        % Generate mpp given the trials
        [mpp_temp] = draw_samples(...
            trials_locations, trials_powers, pi_target_selected, background_rate,...
            v_th_known_related, v_reset_known_related, g_related, gain_related,gamma_related,...
            current_template, funcs, delay_params,stim_threshold,time_max);
        mpp_undefined{iter}=mpp_temp;trials_locations_undefined{iter}=trials_locations;trials_powers_undefined{iter}=trials_powers;
        
        for i_trial = 1:size(cells_probabilities_undefined,1)
            outputs_undefined(i_trial,1)=length(mpp_undefined{iter}(i_trial).times);
        end
        n_trials=n_trials+i_trial;
    end
    %-------
    % Conduct trials on group B, the potentially disconnected cells
    mpp_disconnected{iter}=[];trials_locations_disconnected{iter}=[];trials_powers_disconnected{iter}=[];
    if sum(potentially_disconnected_cells{iter})>0
        % Find cells with close to zero gammas
        cell_list= find(potentially_disconnected_cells{iter});
        gamma_estimates_confirm = 0.5*ones(length(cell_list),1);% for drawing samples...
        [trials_locations,  trials_powers] = random_design(target_locations_selected,power_selected,...
            inner_normalized_products,single_spot_threshold,gamma_estimates_confirm,0,...
            id_notconnected, loc_to_cell,cell_list,n_spots_per_trial,K_disconnected,n_replicates);
        [cells_probabilities_disconnected, ~] = get_prob_and_size(...
            pi_target_selected,trials_locations,trials_powers,stim_unique,prob_trace);
        % Generate mpp given the trials
        [mpp_temp] = draw_samples(...
            trials_locations, trials_powers, pi_target_selected, background_rate,...
            v_th_known_related, v_reset_known_related, g_related, gain_related,gamma_related,...
            current_template,  funcs,    delay_params,stim_threshold,time_max);
        mpp_disconnected{iter}=mpp_temp;trials_locations_disconnected{iter}=trials_locations;trials_powers_disconnected{iter}=trials_powers;
        for i_trial = 1:size(cells_probabilities_disconnected,1)
            outputs_disconnected(i_trial,1)=length(mpp_disconnected{iter}(i_trial).times);
        end
        n_trials=n_trials+i_trial;
    end
    
    %-------
    % Conduct trials on group C, the potentially connected cells
    mpp_connected{iter}=[];
    trials_locations_connected{iter}=[];
    trials_powers_connected{iter}=[];
    if sum(potentially_connected_cells{iter})>0
        % Find cells with close to zero gammas
        cell_list= find(potentially_connected_cells{iter});
        gamma_estimates_confirm = 0.5*ones(length(cell_list),1);% for drawing samples...
        [trials_locations,  trials_powers] = random_design(target_locations_nuclei,power_nuclei,...
            inner_normalized_products,single_spot_threshold,gamma_estimates_confirm,0,...
            connected,  loc_to_cell_nuclei, cell_list,1,K_connected,n_replicates);
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
        cell_list= find(sum(cells_probabilities_undefined>prob_lower_bound,1)>0);
        % Update variational and prior distribution
        variational_params=struct([]);
        for i_cell_idx = 1:length(cell_list)
            i_cell=cell_list(i_cell_idx);
            variational_params(i_cell_idx).pi = variational_params_path.pi(i_cell,iter);
            variational_params(i_cell_idx).p_logit = log(variational_params(i_cell_idx).pi/(1-variational_params(i_cell_idx).pi));
            variational_params(i_cell_idx).alpha = variational_params_path.alpha(i_cell,iter);
            variational_params(i_cell_idx).beta = variational_params_path.beta(i_cell,iter);
        end
        prior_params.pi0= [variational_params(:).pi]';
        prior_params.alpha0= [variational_params(:).alpha]';
        prior_params.beta0 = [variational_params(:).beta]';
        
        designs_remained=cells_probabilities_undefined(:,cell_list);
        active_trials=find(sum(designs_remained,2)>1e-3);
        designs_remained=designs_remained(active_trials,:);
        outputs_remained=outputs_undefined(active_trials,:);
        
        % find neighbours that are not in cell_list:
        % We will treat the gammas in this list as fixed. 
%         neighbour_list=find(sum(cell_neighbours(cell_list,:),1)>0)';
%         neighbour_list=setdiff(neighbour_list,cell_list);
%         designs_neighbours=cells_probabilities_undefined(active_trials,neighbour_list);
%         gamma_neighbours=mean_gamma_current(neighbour_list);
        gamma_neighbours=[];designs_neighbours=[];

        lklh_func=@calculate_likelihood_bernoulli;
        % calculate_likelihood_bernoulli for multiple events
        [parameter_history,~] = fit_working_model_vi(...
            designs_remained,outputs_remained,background_rt, ...
            variational_params,prior_params,C_threshold,...
            designs_neighbours,gamma_neighbours,...
            S,epsilon,eta_logit,eta_beta,maxit,lklh_func);
        
        % Record the variational parameters
        variational_params_path.pi(cell_list,iter+1) = parameter_history.pi(:,end);
        variational_params_path.alpha(cell_list,iter+1) = parameter_history.alpha(:,end);
        variational_params_path.beta(cell_list,iter+1) = parameter_history.beta(:,end);
        
        [mean_gamma_temp, ~] = calculate_posterior_mean(parameter_history.alpha(:,end),parameter_history.beta(:,end),0,1);
        
        mean_gamma_undefined(cell_list,1)=mean_gamma_temp;
        mean_gamma_current(cell_list)=mean_gamma_temp;
        gamma_path(cell_list,iter+1)=mean_gamma_temp;
    end
    %-------------------------------------------------------------%
    
    %----------------------------------------------------------------%
    % Fit the VI on Group B: potentially disconnected cells
    mean_gamma_disconnected=ones(length(related_cell_list),1);
    
    if sum(potentially_disconnected_cells{iter})>0
%         cell_list= find(potentially_disconnected_cells{iter});
        cell_list= find(sum(cells_probabilities_disconnected>prob_lower_bound,1)>0);
        
        variational_params=struct([]);
        for i_cell_idx = 1:length(cell_list)
            i_cell=cell_list(i_cell_idx);
            variational_params(i_cell_idx).pi = variational_params_path.pi(i_cell,iter);
            variational_params(i_cell_idx).p_logit = log(variational_params(i_cell_idx).pi/(1-variational_params(i_cell_idx).pi));
            variational_params(i_cell_idx).alpha = variational_params_path.alpha(i_cell,iter);
            variational_params(i_cell_idx).beta = variational_params_path.beta(i_cell,iter);
        end
        
        prior_params.pi0= [variational_params(:).pi]';
        prior_params.alpha0= [variational_params(:).alpha]';
        prior_params.beta0 = [variational_params(:).beta]';
        
        % Include only the remaining cells
        designs_remained=cells_probabilities_disconnected(:,cell_list);
        active_trials=find(sum(designs_remained,2)>1e-3);
        designs_remained=designs_remained(active_trials,:);
        outputs_remained=outputs_disconnected(active_trials,:);
        
        % find neighbours that are not in cell_list:
%         neighbour_list=find(sum(cell_neighbours(cell_list,:),1)>0)';
%         neighbour_list=setdiff(neighbour_list,cell_list);
%         designs_neighbours=cells_probabilities_disconnected(active_trials,neighbour_list);
%         gamma_neighbours=mean_gamma_current(neighbour_list);
         gamma_neighbours=[];designs_neighbours=[];
         
        lklh_func=@calculate_likelihood_bernoulli;
        [parameter_history,~] = fit_working_model_vi(...
            designs_remained,outputs_remained,background_rt, ...
            variational_params,prior_params,C_threshold,...
            designs_neighbours,gamma_neighbours,...
            S,epsilon,eta_logit,eta_beta,maxit,lklh_func);
        
        variational_params_path.pi(cell_list,iter+1) = parameter_history.pi(:,end);
        variational_params_path.alpha(cell_list,iter+1) = parameter_history.alpha(:,end);
        variational_params_path.beta(cell_list,iter+1) = parameter_history.beta(:,end);
        
        [mean_gamma_temp, ~] = calculate_posterior_mean(parameter_history.alpha(:,end),parameter_history.beta(:,end),0,1);
        
        mean_gamma_disconnected(cell_list,1)=mean_gamma_temp;
        mean_gamma_current(cell_list)=mean_gamma_temp;
        gamma_path(cell_list,iter+1)=mean_gamma_temp;
    end
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
                variational_params(i_cell_idx).pi = variational_params_path.pi(i_cell,iter);
                variational_params(i_cell_idx).p_logit = log(variational_params(i_cell_idx).pi/(1-variational_params(i_cell_idx).pi));
                variational_params(i_cell_idx).alpha = variational_params_path.alpha(i_cell,iter);
                variational_params(i_cell_idx).beta = variational_params_path.beta(i_cell,iter);
                variational_params(i_cell_idx).alpha_gain = variational_params_path.alpha_gain(i_cell,iter);
                variational_params(i_cell_idx).beta_gain = variational_params_path.beta_gain(i_cell,iter);
            end
            
            prior_params.pi0= [variational_params(:).pi]';
            prior_params.alpha0= [variational_params(:).alpha]';
            prior_params.beta0 = [variational_params(:).beta]';
            prior_params.alpha0_gain= [variational_params(:).alpha_gain]';
            prior_params.beta0_gain =[variational_params(:).beta_gain]';
            
            designs_remained=stim_size_connected(active_trials,neighbour_list);
%             active_trials=find(sum(designs_remained,2)>stim_threshold);
%             designs_remained=designs_remained(active_trials,:);
            mpp_remained=mpp_connected{iter}(active_trials);
            
            %             % find neighbours that are not in cell_list:
            %             neighbour_list=find(sum(cell_neighbours(cell_list(cluster_of_cells{i_cluster}),:),1)>0)';
            %             neighbour_list=setdiff(neighbour_list,cell_list(cluster_of_cells{i_cluster}));
            %             designs_neighbours=stim_size_connected(active_trials,neighbour_list);
            %             gamma_neighbours=mean_gamma_current(neighbour_list);
            %               gain_neighbours=mean_gain_current(neighbour_list);
            %
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
            variational_params_path.pi(neighbour_list,iter+1) = parameter_history.pi(:,end);
        variational_params_path.alpha(neighbour_list,iter+1) = parameter_history.alpha(:,end);
        variational_params_path.beta(neighbour_list,iter+1) = parameter_history.beta(:,end);
         variational_params_path.alpha_gain(neighbour_list,iter+1) = parameter_history.alpha_gain(:,end);
        variational_params_path.beta_gain(neighbour_list,iter+1) = parameter_history.beta_gain(:,end);
       
        [mean_gamma_temp, var_gamma_temp] = calculate_posterior_mean(...
            parameter_history.alpha(:,end),parameter_history.beta(:,end),0,1);
        [mean_gain_temp, ~] = calculate_posterior_mean(...
            parameter_history.alpha_gain(:,end),parameter_history.beta_gain(:,end),gain_bound.low,gain_bound.up);
        
            variance_gamma_connected(neighbour_list)=var_gamma_temp;
            var_gamma_path(neighbour_list,iter+1)=var_gamma_temp;
            mean_gamma_connected(neighbour_list,1)=mean_gamma_temp;
            mean_gamma_current(neighbour_list)=mean_gamma_temp;
            mean_gain_current(neighbour_list)=mean_gain_temp;
            gamma_path(neighbour_list,iter+1)=mean_gamma_temp;
        
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
    undefined_to_disconnected=find(sum(cell_neighbours(undefined_to_disconnected,:),1)>0)';
    undefined_to_connected =find(sum(cell_neighbours(undefined_to_connected,:),1)>0);
    % if there are conflicts, move them to the potentially connected cells
    undefined_to_disconnected=setdiff(undefined_to_disconnected,undefined_to_connected);
    
    disconnected_to_undefined = intersect(find(mean_gamma_disconnected>disconnected_confirm_threshold),...
        find(potentially_disconnected_cells{iter}));
    disconnected_to_dead = intersect(find(mean_gamma_disconnected<disconnected_confirm_threshold),...
        find(potentially_disconnected_cells{iter}));
    
    disconnected_to_undefined=find(sum(cell_neighbours(disconnected_to_undefined,:),1)>0);
    % if there are conflicts, move them to the potentially connected cells
    disconnected_to_dead=setdiff(disconnected_to_dead,disconnected_to_undefined);
    
    
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
end
    %%
      [find(gamma_truth(target_cell_list.primary)>0) find(alive_cells{iter})]
  
%% ------------------------------------------ 
    % Plotting 
    
    % Column I: preprocessing
    
    % I.a a 3D plot of cells
    x=cell_locations(:,1);
    y=cell_locations(:,2);
    z=cell_locations(:,3);
    

    figure(101)
    figure(101)
   ax1=axes('Position',[0 0 1 1],'Visible','off');
   ax2=axes('Position',[0.1 0.2 .8 .8],'Visible','off');
 axes(ax2)
    scatter3(-y,x,z,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor','b',...
        'MarkerFaceAlpha',0.6)
    view(-30,10)
    xlim([-150 150]);
    ylim([-150 150]);
    zlim([0 150]);
    axes(ax1)
      title_text = {'I.a all cells'};
      text(0.4,0.08,title_text,'fontsize',15)
 
    saveas(101,strcat('./Figures/work_flow/','FigureIa','.jpg'));

    %%
    % I.b the 3D plot with a plane 
    x=cell_locations(:,1);
    y=cell_locations(:,2);
    z=cell_locations(:,3);
    
    figure(102)
   ax1=axes('Position',[0 0 1 1],'Visible','off');
   ax2=axes('Position',[0.1 0.2 .8 .8],'Visible','off');
   
   sc2=scatter3(ax2,-y,x,z,...
       'MarkerEdgeColor','b',...
       'MarkerFaceColor','b',...
       'MarkerFaceAlpha',0.2);
   hold on;
   sc1=scatter3(ax2,-y(cell_group_list{this_plane}),x(cell_group_list{this_plane}),...
       z(cell_group_list{this_plane}),...
       'MarkerEdgeColor','k',...
       'MarkerFaceColor','k',...
       'MarkerFaceAlpha',0.8);
   hold on;
   view(-30,10)
   xlim([-150 150]);
   ylim([-150 150]);
   zlim([0 150]);
   
   [x_coord, y_coord] = meshgrid(-150:10:150);
   z_coord= mean(cell_locations(cell_group_list{this_plane},3))*ones(size(x_coord,1),size(x_coord,2));
   surf(ax2,-y_coord,x_coord,z_coord,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','none');
   hold on;
   
   title_text = {'I.b selected cells'};
   axes(ax1)
%    text(ax2,-50,0,-50,title_text)
     text(ax2,-40,80,-50,title_text,'fontsize',15)
 
   hold on;
   %scatter3(-50,-50,-50,'MarkerSize',20,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerFaceAlpha',0.8);
   
    hold off;
     saveas(102,strcat('./Figures/work_flow/','FigureIb','.jpg'));

    %%
    % I.c projection of cells onto that plane 
    figure(103)
   ax1=axes('Position',[0 0 1 1],'Visible','off');
   ax2=axes('Position',[0.1 0.2 .8 .8],'Visible','off');
    axes(ax2)
    scatter(-cell_locations(cell_group_list{this_plane},2),...
        cell_locations(cell_group_list{this_plane},1),...
        'Marker','o','SizeData',60,...
        'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.5)
    hold on;
    scatter(-target_locations_selected(:,2),...
        target_locations_selected(:,1),...
        'Marker','d','SizeData',60,...
        'MarkerFaceColor','r', 'MarkerEdgeColor','r', 'MarkerFaceAlpha',0.5)
    
     xlim([-160 160]);
    ylim([-160 160]);
    axis off;
    
    axes(ax1)
    title_text = {'I.c projected cells and selected locations'};
    axes(ax1)
    text(0.15,0.08,title_text,'fontsize',15)
   
    hold off;
    
    saveas(103,strcat('./Figures/work_flow/','FigureIc','.jpg'));

    
    %% ---------------------------------------
    % Column II: heatmaps of trials
    iter_II=4;
    find(potentially_disconnected_cells{iter_II})
    find(potentially_connected_cells{iter_II})
    find(undefined_cells{iter_II})
    %%
    % II.a trials for undeifned cells 
    figure(201)
    ax1=axes('Position',[0 0 1 1],'Visible','off');
   ax2=axes('Position',[0.1 0.2 .8 .8],'Visible','off');
    axes(ax2)
    
           for i_cell = 1:length(cell_group_list{this_plane})
                cell_index =cell_group_list{this_plane}(i_cell);
                
                    if undefined_cells{iter_II}(i_cell)==1
                        facecolor='k';markershape='o';markeralpha=0.4;
                    elseif dead_cells{iter_II}(i_cell)==1
                        facecolor='r';markershape='x';markeralpha=0.2;
                    elseif alive_cells{iter_II}(i_cell)==1
                      facecolor='b';markershape='o';markeralpha=0.2;
                    elseif potentially_connected_cells{iter_II}(i_cell)==1
                      facecolor='g';markershape='o';markeralpha=0.2;
                    elseif potentially_disconnected_cells{iter_II}(i_cell)==1
                      facecolor='r';markershape='o';markeralpha=0.2;
                    end
                    
                    scatter(cell_locations(cell_index,2),...
                        cell_locations(cell_index,1),...
                        'Marker',markershape,'SizeData',150,...
                        'MarkerFaceColor',facecolor, 'MarkerEdgeColor',facecolor,...
                        'MarkerFaceAlpha',markeralpha)
                  
                
                hold on;
            end
            %xlabel(strcat(num2str(n_total_trial(iter)), ' trials'));
            xlim([-160 160]);
ylim([-160 160]);

facecolor='k';
    for i_trial = 1:size(trials_locations_undefined{iter_II},1)
           scatter(target_locations_selected(trials_locations_undefined{iter_II}(i_trial,:),2),...
                        target_locations_selected(trials_locations_undefined{iter_II}(i_trial,:),1),...
                        'Marker','d','SizeData',60,...
                        'MarkerFaceColor',facecolor, 'MarkerEdgeColor',facecolor,...
                        'MarkerFaceAlpha',0.2)
    end
     axis off;
    
    axes(ax1)
    title_text = {'II.a designed trials on undefined cells'};
    axes(ax1)
   text(0.2,0.08,title_text,'fontsize',15)
   
%      target_locations_selected
 saveas(201,strcat('./Figures/work_flow/','FigureIIa','.jpg'));

    
    %%
    % II.b trials for disconnected cells
    figure(202)
       ax1=axes('Position',[0 0 1 1],'Visible','off');
   ax2=axes('Position',[0.1 0.2 .8 .8],'Visible','off');
    axes(ax2)
    for i_cell = 1:length(cell_group_list{this_plane})
        cell_index =cell_group_list{this_plane}(i_cell);
        
        if undefined_cells{iter_II}(i_cell)==1
            facecolor='k';markershape='o';markeralpha=0.2;
        elseif dead_cells{iter_II}(i_cell)==1
            facecolor='r';markershape='x';markeralpha=0.2;
        elseif alive_cells{iter_II}(i_cell)==1
            facecolor='b';markershape='o';markeralpha=0.2;
        elseif potentially_connected_cells{iter_II}(i_cell)==1
            facecolor='g';markershape='o';markeralpha=0.2;
        elseif potentially_disconnected_cells{iter_II}(i_cell)==1
            facecolor='r';markershape='o';markeralpha=0.4;
        end
        
        scatter(cell_locations(cell_index,2),...
            cell_locations(cell_index,1),...
            'Marker',markershape,'SizeData',150,...
            'MarkerFaceColor',facecolor, 'MarkerEdgeColor',facecolor,...
            'MarkerFaceAlpha',markeralpha)
        
        
        hold on;
    end
    %xlabel(strcat(num2str(n_total_trial(iter)), ' trials'));
    xlim([-160 160]);
    ylim([-160 160]);
    axis off;
    facecolor='r';
    for i_trial = 1:size(trials_locations_disconnected{iter_II},1)
        scatter(target_locations_selected(trials_locations_disconnected{iter_II}(i_trial,:),2),...
            target_locations_selected(trials_locations_disconnected{iter_II}(i_trial,:),1),...
            'Marker','d','SizeData',60,...
            'MarkerFaceColor',facecolor, 'MarkerEdgeColor',facecolor,...
            'MarkerFaceAlpha',0.2)
    end
    %      target_locations_selected
      axes(ax1)
    title_text = {'II.b designed trials on pot. disconnected cells'};
    axes(ax1)
   % text(0.25,0.05,title_text)
  text(0.2,0.08,title_text,'fontsize',15)
   
     saveas(202,strcat('./Figures/work_flow/','FigureIIb','.jpg'));


    %%
    % II.c trials for connected cells (3D)
 figure(203)
    ax1=axes('Position',[0 0 1 1],'Visible','off');
   ax2=axes('Position',[0.1 0.2 .8 .8],'Visible','off');
    axes(ax2)
    for i_cell = 1:length(cell_group_list{this_plane})
        cell_index =cell_group_list{this_plane}(i_cell);
        
        if undefined_cells{iter_II}(i_cell)==1
            facecolor='k';markershape='o';markeralpha=0.2;
        elseif dead_cells{iter_II}(i_cell)==1
            facecolor='r';markershape='x';markeralpha=0.2;
        elseif alive_cells{iter_II}(i_cell)==1
            facecolor='b';markershape='o';markeralpha=0.2;
        elseif potentially_connected_cells{iter_II}(i_cell)==1
            facecolor='g';markershape='o';markeralpha=0.5;
        elseif potentially_disconnected_cells{iter_II}(i_cell)==1
            facecolor='r';markershape='o';markeralpha=0.2;
        end
        
        scatter3(-cell_locations(cell_index,2),...
            cell_locations(cell_index,1),cell_locations(cell_index,3),...
            'Marker',markershape,'SizeData',50,...
            'MarkerFaceColor',facecolor, 'MarkerEdgeColor',facecolor,...
            'MarkerFaceAlpha',markeralpha)
        hold on;
    end
    facecolor='g';
    for i_trial = 1:size(trials_locations_connected{iter_II},1)
           scatter3(-target_locations_nuclei(trials_locations_connected{iter_II}(i_trial,:),2),...
                        target_locations_nuclei(trials_locations_connected{iter_II}(i_trial,:),1),...
                        target_locations_nuclei(trials_locations_connected{iter_II}(i_trial,:),3),...
                        'Marker','d','SizeData',50,...
                        'MarkerFaceColor',facecolor, 'MarkerEdgeColor',facecolor,...
                        'MarkerFaceAlpha',0.2)
    end
    %xlabel(strcat(num2str(n_total_trial(iter)), ' trials'));
    view(-30,40)
    xlim([-150 150]);
    ylim([-150 150]);
    zlim([45 65]);
      title_text = {'II.c designed trials on pot. connected cells'};
   axes(ax1)
   text(ax2,-150,40,30,title_text,'fontsize',15)
   %text(0.15,0.08,title_text,'fontsize',15)
   
   hold on;
 
     saveas(203,strcat('./Figures/work_flow/','FigureIIc','.jpg'));

    %% ---------------------------------------
    % Column III: observed data 
    %%
    % III.a Response and Probability 
    figure(301)
        ax1=axes('Position',[0 0 1 1],'Visible','off');
   ax2=axes('Position',[0.1 0.2 .8 .8],'Visible','off');
    axes(ax2)
    [cells_probabilities_undefined, ~] = get_prob_and_size(...
        pi_target_selected,trials_locations_undefined{iter_II},trials_powers_undefined{iter_II},...
        stim_unique,prob_trace);
    
    for i_trial = 1:size(cells_probabilities_undefined,1)
        outputs_undefined(i_trial,1)=1*(1-isempty(mpp_undefined{iter_II}(i_trial).times));
    end
    
        cell_list=find(sum(cells_probabilities_undefined,1)>0.1);
    facecolor='k';
    for i_trial = 1:size(trials_locations_undefined{iter_II},1)
        scatter(-1,...
            i_trial,...
            'Marker','s','SizeData',20,...
            'MarkerFaceColor','b', 'MarkerEdgeColor','w',...
            'MarkerFaceAlpha',max(outputs_undefined(i_trial),0.01))
        hold on;
       for i_cell = 1:length(cell_list)
           i_cell_idx=cell_list(i_cell);
           if cells_probabilities_undefined(i_trial,i_cell_idx)>0.01
       scatter(i_cell+4,...
            i_trial,...
            'Marker','s','SizeData',20,...
            'MarkerFaceColor',facecolor, 'MarkerEdgeColor','w',...
            'MarkerFaceAlpha',min(max(cells_probabilities_undefined(i_trial,i_cell_idx),0.01),1))
           end
       end
    end
    %      target_locations_selected
      xlim([-3 i_cell+5]);
      xticks([-1 3  4+(1:length(cell_list)) ])
      xticklabels({'Response' 'Cells:'  string(cell_list)})
    ylim([0 i_trial]);
    ylabel('Trials');
  hold on;
      title_text = {'III.a responses and firing probabilities (undefined)'};
   axes(ax1)
   text(0.15,0.08,title_text,'fontsize',15)
   hold off;
 
   saveas(301,strcat('./Figures/work_flow/','FigureIIIa','.jpg'));

  % Fitted values
%   cell_list=find(undefined_cells{iter_II});
%   alpha_current=exp(variational_params_path.log_alpha(cell_list,iter_II));
%   beta_current=exp(variational_params_path.log_beta(cell_list,iter_II));
%   figure(3012)
%  xgrid=0:0.02:1; 
%     facecolor='k';
%        for i_cell = 1:length(cell_list)
%            
%        line(xgrid,betapdf(xgrid,alpha_current(i_cell),beta_current(i_cell)));
%        end
%     
%     %      target_locations_selected
%       xlim([0 1]);
%     ylim([0 5]);
%   hold off;
%        saveas(3012,strcat('./Figures/work_flow/','FigureIIIa2','.jpg'));
%%
  % III.b Response and Probability (disconnected)
   [cells_probabilities_disconnected, ~] = get_prob_and_size(...
        pi_target_selected,trials_locations_disconnected{iter_II},trials_powers_disconnected{iter_II},...
        stim_unique,prob_trace);
    outputs_disconnected(:)=[];
    for i_trial = 1:size(cells_probabilities_disconnected,1)
        outputs_disconnected(i_trial,1)=1*(1-isempty(mpp_disconnected{iter_II}(i_trial).times));
    end
    
  figure(302)
    ax1=axes('Position',[0 0 1 1],'Visible','off');
   ax2=axes('Position',[0.1 0.2 .8 .8],'Visible','off');
    axes(ax2)
   
    cell_list=find(sum(cells_probabilities_disconnected,1)>0.1);
    facecolor='k';
    for i_trial = 1:size(trials_locations_disconnected{iter_II},1)
        scatter(-1,...
            i_trial,...
            'Marker','s','SizeData',20,...
            'MarkerFaceColor','b', 'MarkerEdgeColor','w',...
            'MarkerFaceAlpha',max(outputs_disconnected(i_trial),0.01))
        hold on;
           
       for i_cell = 1:length(cell_list)
            i_cell_idx=cell_list(i_cell);
           if cells_probabilities_disconnected(i_trial,i_cell_idx)>0.01
       scatter(1.3*i_cell+4,...
            i_trial,...
            'Marker','s','SizeData',20,...
            'MarkerFaceColor',facecolor, 'MarkerEdgeColor','w',...
            'MarkerFaceAlpha',min(max(cells_probabilities_disconnected(i_trial,i_cell_idx),0.01),1))
           end
       end
    end
    %      target_locations_selected
    xlim([-3 1.3*i_cell+5]);
      xticks([-1 3  4+1.3*(1:length(cell_list)) ])
      xticklabels({'Response' 'Cells:'  string(cell_list)})
        ylim([0 i_trial]);
    ylabel('Trials');
  hold on;
 
   title_text = {'III.b responses and firing probabilities (disconnected)'};
   axes(ax1)
   text(0.1,0.08,title_text,'fontsize',15)
 
  hold off;
  saveas(302,strcat('./Figures/work_flow/','FigureIIIb','.jpg'));

    %%
    % III.b Process and Stim size 
    
 figure(303)
 ax1=axes('Position',[0 0 1 1],'Visible','off');
 ax2=axes('Position',[0.1 0.2 .8 .8],'Visible','off');
 axes(ax2)
 
      [~, stim_size_connected] = get_prob_and_size(...
            pi_target_nuclei,trials_locations_connected{iter_II},trials_powers_connected{iter_II},...
            stim_unique,prob_trace);
           
        transparencies= stim_size_connected/max(max(stim_size_connected));
        cell_list=find(sum(transparencies)>4);
    facecolor='k';
    
    for i_trial = 1:size(trials_locations_connected{iter_II},1)
        if ~isempty(mpp_connected{iter_II}(i_trial).times)
           for i_event = 1:length(mpp_connected{iter_II}(i_trial).times)
               x_sp=[mpp_connected{iter_II}(i_trial).times(i_event) mpp_connected{iter_II}(i_trial).times(i_event)];
               x_sp=x_sp/(time_max/10);
               y_sp=[i_trial+0.2 i_trial+0.8];
               line(x_sp,y_sp,'LineWidth',2)
               
           end
        end
          
        
       for i_cell = 1:length(cell_list)
           i_cell_idx = cell_list(i_cell);
           %if transparencies(i_trial,i_cell)>0.2
       scatter(i_cell+11,...
            i_trial,...
            'Marker','s','SizeData',20,...
            'MarkerFaceColor',facecolor, 'MarkerEdgeColor','w',...
            'MarkerFaceAlpha', transparencies(i_trial,i_cell_idx))
       
        hold on;
           %end
       end
    end
       xlim([-2 i_cell+12]);
      xticks([1 4 6 8 10.5  11+(1:length(cell_list)) ])
      xticklabels({'Time (ms):' string([3 9 12]) 'Cells:'  string(cell_list)})
        ylim([0 i_trial]);
    ylabel('Trials');
  hold on;
 
   title_text = {'III.c responses and firing probabilities (connected)'};
   axes(ax1)
   
   text(0.15,0.08,title_text,'fontsize',15)
 
    gamma_truth(find(potentially_connected_cells{iter_II}));
       saveas(303,strcat('./Figures/work_flow/','FigureIIIc','.jpg'));

    %% Draw the estimated gamma and gains 
%       figure(3022)
%   % Fitted values
%   cell_list=find(potentially_connected_cells{iter_II});
%   alpha_current=exp(variational_params_path.log_alpha(cell_list,iter_II));
%   beta_current=exp(variational_params_path.log_beta(cell_list,iter_II));
%   
%  xgrid=0:0.02:1; 
%   figure(4)
%     facecolor='k';
%        for i_cell = 1:length(cell_list)
%            
%        line(xgrid,betapdf(xgrid,alpha_current(i_cell),beta_current(i_cell)));
%        end
%     
%     %      target_locations_selected
%       xlim([0 1]);
%     ylim([0 5]);
%   hold off;
%        saveas(3022,strcat('./Figures/work_flow/','FigureIIIb2','.jpg'));

  %% Gains 
%   alpha_gain_current=exp(variational_params_path.log_alpha(cell_list,iter_II));
%   beta_gain_current=exp(variational_params_path.log_beta(cell_list,iter_II));
%   
%   figure(3023)
%   xgrid=0:0.02:1;
%   scaled_xgrid=(gain_bound.low+ (gain_bound.up-gain_bound.low)*xgrid);
%   figure(5)
%   facecolor='k';
%   for i_cell = 1:length(cell_list)
%       
%       line(scaled_xgrid,betapdf(xgrid,alpha_current(i_cell),beta_current(i_cell)));
%       
%   end
%   
%   %      target_locations_selected
%   xlim([0.01 0.03]);
%   ylim([0 5]);
%   hold off;
%      saveas(3023,strcat('./Figures/work_flow/','FigureIIIb3','.jpg'));

           
  
    %% ---------------------------------------
    % Column IV: assignments
    % IV.a undefined  cells  
    
    
    cell_list=find(undefined_cells{iter_II});
    if vf_type == 1
        alpha_current=exp(variational_params_path.log_alpha(cell_list,iter_II));
        beta_current=exp(variational_params_path.log_beta(cell_list,iter_II));
        gamma_mean = alpha_current./(alpha_current+beta_current);
        
    elseif vf_type==2
        alpha_current=variational_params_path.log_alpha(cell_list,iter_II);
        beta_current=exp(variational_params_path.log_beta(cell_list,iter_II));
        gamma_mean= zeros(length(cell_list),1);gamma_var= zeros(length(cell_list),1);
        for i_temp = 1:length(cell_list)
            gamma_i=beta_current(i_temp)*normal_samples+alpha_current(i_temp);
            gamma_i=exp(gamma_i)./(1+exp(gamma_i));
            gamma_mean(i_temp)=mean(gamma_i);
            gamma_var(i_temp)=var(gamma_i);
        end
        
    end
    xgrid=0:0.02:1;
    
    figure(401)
    ax1=axes('Position',[0 0 1 1],'Visible','off');
    ax2=axes('Position',[0.1 0.2 .8 .8],'Visible','off');
    axes(ax2)
    for i_cell = 1:length(cell_list)
        if gamma_mean(i_cell)>connected_threshold
            color_codes = [0 1 0 0.4]; %green
        elseif gamma_mean(i_cell)<disconnected_threshold
            color_codes = [1 0 0 0.4]; %red
        else
            color_codes = [0 0 0 0.2]; %black 
        end
        line([gamma_mean(i_cell) gamma_mean(i_cell)],[4 5],'LineStyle','-','LineWidth',1,'Color',color_codes)
        if vf_type == 1
        line(xgrid,betapdf(xgrid,alpha_current(i_cell),beta_current(i_cell)),'Color',color_codes );
        elseif vf_type == 2
        xlogit=log(xgrid./(1-xgrid));
            line(xgrid,normpdf(xlogit,alpha_current(i_cell),beta_current(i_cell)),'Color',color_codes );
            
        end 
        hold on;
    end
    line([connected_threshold connected_threshold],[3 5],'LineStyle','-','LineWidth',1,'Color',[0 1 0])
    line([disconnected_threshold disconnected_threshold],[3 5],'LineStyle','-','LineWidth',1,'Color',[1 0 0])
    
    %      target_locations_selected
    xlim([0 1]);
    ylim([0 5]);
    axis on;
    ylabel('Density');
    xlabel('\gamma and E(\gamma)')
    hold on;
    
    title_text = {'IV.a posterior dist. of \gamma (undefined)'};
    axes(ax1)
   text(0.23,0.06,title_text,'fontsize',15)
    hold off;
       saveas(401,strcat('./Figures/work_flow/','FigureIVa','.jpg'));

       %%
         % IV.b potentially disconnected cells  
    cell_list=find(potentially_disconnected_cells{iter_II});
    
     if vf_type == 1
        alpha_current=exp(variational_params_path.log_alpha(cell_list,iter_II));
        beta_current=exp(variational_params_path.log_beta(cell_list,iter_II));
        gamma_mean = alpha_current./(alpha_current+beta_current);     
    elseif vf_type==2
        alpha_current=variational_params_path.log_alpha(cell_list,iter_II);
        beta_current=exp(variational_params_path.log_beta(cell_list,iter_II));
        gamma_mean= zeros(length(cell_list),1);gamma_var= zeros(length(cell_list),1);
        for i_temp = 1:length(cell_list)
            gamma_i=beta_current(i_temp)*normal_samples+alpha_current(i_temp);
            gamma_i=exp(gamma_i)./(1+exp(gamma_i));
            gamma_mean(i_temp)=mean(gamma_i);
            gamma_var(i_temp)=var(gamma_i);
        end
        
    end
    xgrid=0:0.02:1;
    figure(402)
    ax1=axes('Position',[0 0 1 1],'Visible','off');
    ax2=axes('Position',[0.1 0.2 .8 .8],'Visible','off');
    axes(ax2)
    for i_cell = 1:length(cell_list)
        if gamma_mean(i_cell)<disconnected_threshold
            color_codes = [1 0 0 0.4]; %red
        else
            color_codes = [0 0 0 0.2]; %black 
        end
        line([gamma_mean(i_cell) gamma_mean(i_cell)],[4 5],'LineStyle','-','LineWidth',1,'Color',color_codes)
         if vf_type == 1
        line(xgrid,betapdf(xgrid,alpha_current(i_cell),beta_current(i_cell)),'Color',color_codes );
        elseif vf_type == 2
        xlogit=log(xgrid./(1-xgrid));
            line(xgrid,normpdf(xlogit,alpha_current(i_cell),beta_current(i_cell)),'Color',color_codes );
            
        end 
        hold on;
    end
    line([disconnected_threshold disconnected_threshold],[3 5],'LineStyle','-','LineWidth',1,'Color',[1 0 0])
    
    %      target_locations_selected
     xlim([0 1]);
    ylim([0 5]);
    axis on;
    ylabel('Density');
    xlabel('\gamma and E(\gamma)')
    hold on;
    
    title_text = {'IV.b posterior dist. of \gamma (disconnected)'};
    axes(ax1)
    
   text(0.2,0.06,title_text,'fontsize',15)
 
       saveas(402,strcat('./Figures/work_flow/','FigureIVb','.jpg'));

       
  %%
    % IV.b Potentially connected cells 
      % Fitted values
  cell_list=find(potentially_connected_cells{iter_II});
  
  if vf_type == 1
      alpha_current=exp(variational_params_path.log_alpha(cell_list,iter_II));
      beta_current=exp(variational_params_path.log_beta(cell_list,iter_II));
      gamma_mean = alpha_current./(alpha_current+beta_current);
      
      
      %cell_list_old=find(potentially_connected_cells{iter_II-1});
      alpha_current_old=exp(variational_params_path.log_alpha(cell_list,iter_II-1));
      beta_current_old=exp(variational_params_path.log_beta(cell_list,iter_II-1));
      gamma_mean_old = alpha_current_old./(alpha_current_old+beta_current_old);
      
      gamma_change = abs(gamma_mean-gamma_mean_old);
  elseif vf_type==2
      alpha_current=variational_params_path.log_alpha(cell_list,iter_II);
      beta_current=exp(variational_params_path.log_beta(cell_list,iter_II));
      gamma_mean= zeros(length(cell_list),1);gamma_var= zeros(length(cell_list),1);
      for i_temp = 1:length(cell_list)
          gamma_i=beta_current(i_temp)*normal_samples+alpha_current(i_temp);
          gamma_i=exp(gamma_i)./(1+exp(gamma_i));
          gamma_mean(i_temp)=mean(gamma_i);
          gamma_var(i_temp)=var(gamma_i);
      end
      gamma_change=sqrt(gamma_var);
      
  end
    
  xgrid=0:0.02:1; 
  figure(403)
    ax1=axes('Position',[0 0 1 1],'Visible','off');
    ax2=axes('Position',[0.1 0.2 .8 .8],'Visible','off');
    facecolor='k';
       for i_cell = 1:length(cell_list)
           
        if gamma_mean(i_cell)>connected_threshold & gamma_change(i_cell)<change_threshold
            color_codes = [0 0 1 0.4]; %blue
        elseif gamma_mean(i_cell)<disconnected_threshold
            color_codes = [1 0 0 0.4]; %red
        else
            color_codes = [0 1 0 0.4]; %green
        end
       line([gamma_mean(i_cell) gamma_mean(i_cell)],[4 5],'LineStyle','-','LineWidth',1,'Color',color_codes)
       
        if vf_type == 1
       line(xgrid,betapdf(xgrid,alpha_current(i_cell),beta_current(i_cell)),'Color',color_codes );
       %line(xgrid,betapdf(xgrid,alpha_current_old(i_cell),beta_current_old(i_cell)),...
       %    'LineStyle','--','Color',color_codes );
        elseif vf_type == 2
        xlogit=log(xgrid./(1-xgrid));
            line(xgrid,normpdf(xlogit,alpha_current(i_cell),beta_current(i_cell)),'Color',color_codes );
            
        end 
       hold on;
       end
    line([connected_threshold connected_threshold],[3 5],'LineStyle','-','LineWidth',1,'Color',[0 0 1])
    line([disconnected_threshold disconnected_threshold],[3 5],'LineStyle','-','LineWidth',1,'Color',[1 0 0])
    
    %      target_locations_selected
      xlim([0 1]);
    ylim([0 5]);
    axis on;
    ylabel('Density');
    xlabel('\gamma and E(\gamma)')
    hold on;
    
    title_text = {'IV.c posterior dist. of \gamma (connected)'};
    axes(ax1)
    
   text(0.2,0.06,title_text,'fontsize',15)
       saveas(403,strcat('./Figures/work_flow/','FigureIVc','.jpg'));

       %%
       
       %% Draw the plots:
final_iter=iter;
%% Final fits

  if vf_type == 1
        alpha_current=exp(variational_params_path.log_alpha(:,final_iter));
        beta_current=exp(variational_params_path.log_beta(:,final_iter));
        gamma_mean = alpha_current./(alpha_current+beta_current);     
    elseif vf_type==2
        alpha_current=variational_params_path.log_alpha(:,final_iter);
        beta_current=exp(variational_params_path.log_beta(:,final_iter));
        gamma_mean= zeros(n_cell_this_plane,1);
        for i_temp = 1:n_cell_this_plane
            gamma_i=beta_current(i_temp)*normal_samples+alpha_current(i_temp);
            gamma_i=exp(gamma_i)./(1+exp(gamma_i));
            gamma_mean(i_temp)=mean(gamma_i);
        end
        
        alpha_current=variational_params_path.log_alpha_gain(:,final_iter);
        beta_current=exp(variational_params_path.log_beta_gain(:,final_iter));
        gain_mean= zeros(n_cell_this_plane,1);
        for i_temp = 1:n_cell_this_plane
            gain_i=beta_current(i_temp)*normal_samples+alpha_current(i_temp);
            gain_i=exp(gain_i)./(1+exp(gain_i));
            gain_i = gain_i*(gain_bound.up-gain_bound.low) +gain_bound.low;
            gain_mean(i_temp)=mean(gain_i);
        end
        
        
  end
  gamma_mean(gamma_mean<disconnected_threshold)=0;
  
  %%
    
figure(1)
scatter(gamma_truth,gamma_mean,'Marker','o','SizeData',25,...
    'MarkerFaceColor','b', 'MarkerEdgeColor','b', 'MarkerFaceAlpha',0.8)
x=[0 1];y=[0 1];
hold on;
line(x,y,'Color','red','LineStyle','--')
hold off;

xlim([0 1]);
ylim([0 1]);

xlabel('True synaptic success rate');
ylabel('Estimated synaptic success rate');
saveas(1,strcat('./Figures/Toy/','Gamma_truth_vs_fits','.jpg'));


%%
figure(2)
scatter(gain_truth(find(gamma_truth>0)),gain_mean(find(gamma_truth>0)),'Marker','o','SizeData',25,...
    'MarkerFaceColor','g', 'MarkerEdgeColor','g', 'MarkerFaceAlpha',0.8)
x=[0 1];y=[0 1];
hold on;
line(x,y,'Color','red','LineStyle','--')
hold off;

xlim([0.01 0.03]);
ylim([0.01 0.03]);

xlabel('True optical gain');
ylabel('Estimated optical gain');
saveas(2,strcat('./Figures/Toy/','Gain_truth_vs_fits','.jpg'));

%% Number of trials v.s. dead cells and alive cells
n_dead_cells = zeros(final_iter,1);
n_alive_cells = zeros(final_iter,1);
n_trial_proc=zeros(final_iter,1);
for i=2:final_iter
    n_dead_cells(i)=sum(dead_cells{i});
    n_alive_cells(i)=sum(alive_cells{i});
    n_trial_proc(i)=length(mpp_undefined{i-1})+length(mpp_connected{i-1})+length(mpp_disconnected{i-1});
end
n_total_trial=cumsum(n_trial_proc);

%%
figure(3)

plot(n_total_trial,n_dead_cells,'Color','red','LineStyle','-','LineWidth',3)

hold on;
line([15 80],[50 50],'Color','red','LineStyle','-','LineWidth',3)
text(85,50,'Number of disconnected cells')

line(n_total_trial,n_alive_cells,'Color','blue','LineStyle','-','LineWidth',3)
line([15 80],[43 43],'Color','blue','LineStyle','-','LineWidth',3)
text(85,43,'Number of connected cells')
xlim([0 max(n_total_trial)]);
ylim([0 n_cell_this_plane]);

xlabel('Number of trials');
ylabel('Cell counts');

hold off;

saveas(3,strcat('./Figures/Toy/','Cell_counts_vs_ntrial','.jpg'));


%% Plot the change over course
final_iter = length(alive_cells);
for iter = 1:final_iter
    figure(iter)
    scatter(cell_locations(cell_group_list{this_plane},2),...
        cell_locations(cell_group_list{this_plane},1),...
        'Marker','o','SizeData',1,...
        'MarkerFaceColor','r', 'MarkerEdgeColor','r', 'MarkerFaceAlpha',0.2)
    hold on;
    for i_cell = 1:length(cell_group_list{this_plane})
        cell_index =cell_group_list{this_plane}(i_cell);
        
        if undefined_cells{iter}(i_cell)==1
            facecolor='k';markershape='o';
        elseif dead_cells{iter}(i_cell)==1
            facecolor='r';markershape='x';
        elseif alive_cells{iter}(i_cell)==1
            facecolor='b';markershape='o';
        elseif potentially_connected_cells{iter}(i_cell)==1
            facecolor='g';markershape='o';
        elseif potentially_disconnected_cells{iter}(i_cell)==1
            facecolor='r';markershape='o';
        end
        
        scatter(cell_locations(cell_index,2),...
            cell_locations(cell_index,1),...
            'Marker',markershape,'SizeData',200,...
            'MarkerFaceColor',facecolor, 'MarkerEdgeColor',facecolor,...
            'MarkerFaceAlpha',0.6)
        
        
        hold on;
    end
    xlabel(strcat(num2str(n_total_trial(iter)), ' trials'));
    xlim([-160 160]);
    ylim([-160 160]);
    
    %xlabel('X (um)');
    %ylabel('Y (um)');
    %             axis off;
    
    
    txt=cell([5 1]);
    txt{1}='Undefined';
    txt{2}='Disconnected';
    txt{3}='Connected';
    txt{4}='Pot. connected';
    txt{5}='Pot. disconnected';
    xcord=[-150 -150 -150 -150 -150];
    ycord=[150 130 110 90 70];
    for i=1:5
        if i==1
            facecolor='k';markershape='o';
        elseif i==2
            facecolor='r';markershape='x';
        elseif i==3
            facecolor='b';markershape='o';
        elseif i==4
            facecolor='g';markershape='o';
        elseif i==5
            facecolor='r';markershape='o';
        end
        scatter(xcord(i),ycord(i),...
            'Marker',markershape,'SizeData',200,...
            'MarkerFaceColor',facecolor, 'MarkerEdgeColor',facecolor,...
            'MarkerFaceAlpha',0.6)
        text(xcord(i)+10,ycord(i),txt{i})
    end
    hold off;
    
%     saveas(iter,strcat('./Figures/Toy/','Cell_map', num2str(iter),'.jpg'));
    
end

%% Figure showing the movement of cells 

% Summarize the experiment:

cell_assignments = zeros(final_iter,n_cell_this_plane);
group_assignments = zeros(final_iter,2,5);
cell_trial_numbers =zeros(final_iter,n_cell_this_plane); % number of trials on each cell
cell_gamma_mean = zeros(final_iter,n_cell_this_plane);
cell_gamma_variance = zeros(final_iter,n_cell_this_plane);
cell_gain_mean = zeros(final_iter,n_cell_this_plane);
cell_gain_variance = zeros(final_iter,n_cell_this_plane);

trial_numbers = zeros(final_iter,3);

connected_ind=find(gamma_truth(cell_group_list{this_plane})>0);
disconnected_ind=find(gamma_truth(cell_group_list{this_plane})==0);

for iter=2:final_iter
    % count the assignments 
    cell_assignments(iter,:)= 5*alive_cells{iter}+4*potentially_connected_cells{iter}+...
        +undefined_cells{iter}*3+potentially_disconnected_cells{iter}*2+1*dead_cells{iter};
   for j = 1:5
    group_assignments(iter,2,j)= sum(cell_assignments(iter,connected_ind)==j);
    group_assignments(iter,1,j)= sum(cell_assignments(iter,disconnected_ind)==j);
    end
    
    % count the number of trials 
    
    % map the stim location to cells 
    
    for j=1:n_cell_this_plane
    cell_trial_numbers(iter,j)= sum(loc_to_cell_nuclei([trials_locations_connected{iter-1}])==j)+...
        sum(sum(trials_locations_undefined{iter-1}==j))+sum(sum(trials_locations_disconnected{iter-1}==j));
    end
    trial_numbers(iter,3)=size(trials_locations_connected{iter-1},1);
    trial_numbers(iter,2)=size(trials_locations_undefined{iter-1},1);
    trial_numbers(iter,1)=size(trials_locations_disconnected{iter-1},1);
    
    % record the gamma and gain 
    for j=1:n_cell_this_plane
        
        [cell_gain_mean(iter,j), cell_gain_variance(iter,j)]=calculate_posterior_mean(...
            variational_params_path.alpha_gain(j,iter),variational_params_path.beta_gain(j,iter),0,1);
        [cell_gamma_mean(iter,j), cell_gamma_variance(iter,j)]=calculate_posterior_mean(...
            variational_params_path.alpha(j,iter),variational_params_path.beta(j,iter),0,1);
     end
end

% total number of trials
total_trial_numbers = cumsum(sum(trial_numbers,2));


%% Draw the change of cell assignments

figure(1)
color_list= [[1 0 0 1]; [1 0 0 0.3]; [0 0 0 0.3]; [0 0 1 1]; [0 1 0 1]];

x_loc=ones(n_cell_this_plane,1)*0;
y_loc= (1:n_cell_this_plane);
scatter(x_loc,y_loc,...
    'Marker','o','SizeData',20,...
    'MarkerFaceColor',color_list(3,1:3), 'MarkerEdgeColor',color_list(3,1:3),...
    'MarkerFaceAlpha',color_list(3,4))
hold on;
size(group_assignments);
for i = 1:size(group_assignments,1)
    n_trial_temp=total_trial_numbers(i);
    cell_counts = 0;
    for j=1:5 % disconnected
        n_cell_temp=group_assignments(i,1,j);
        if n_cell_temp>0
        x_loc=ones(n_cell_temp,1)*n_trial_temp;
        y_loc=cell_counts + (1:n_cell_temp);
        scatter(x_loc,y_loc,...
            'Marker','o','SizeData',20,...
            'MarkerFaceColor',color_list(j,1:3), 'MarkerEdgeColor',color_list(j,1:3),...
           'MarkerFaceAlpha',color_list(j,4))
        cell_counts=cell_counts+n_cell_temp;
        end
    end
     for j=1:5 %connected
        n_cell_temp=group_assignments(i,2,j);
        if n_cell_temp>0
        x_loc=ones(n_cell_temp,1)*n_trial_temp;
        y_loc=cell_counts + (1:n_cell_temp);
        scatter(x_loc,y_loc,...
            'Marker','o','SizeData',20,...
            'MarkerFaceColor',color_list(j,1:3), 'MarkerEdgeColor',color_list(j,1:3),...
           'MarkerFaceAlpha',color_list(j,4))
        cell_counts=cell_counts+n_cell_temp;
        end
    end
end

%% Draw the number of trials per cell 

figure(1)
hold on;
for j=1:size(cell_trial_numbers,2)
    
    n_trial_temp=0;
    for i=1:(size(cell_trial_numbers,1)-1)
        n_trial_this=cell_trial_numbers(i,j);
        line([total_trial_numbers(i) total_trial_numbers(i+1)],...
            [n_trial_temp n_trial_temp+n_trial_this],...
            'Color',color_list(cell_assignments(i+1,j),:),...
            'LineStyle','-','LineWidth',3)
        n_trial_temp=n_trial_temp+n_trial_this;
    end
    
end

%% Draw trials by groups:
figure(1)
hold on;
for j=1:3
    trials_temp = cumsum(trial_numbers(:,j));
  line(total_trial_numbers,...
             trials_temp,...
            'Color',color_list(j+1,:),...
            'LineStyle','-','LineWidth',3)
end

%% Now draw the change of gamma & gain for the disconnected cells 



figure(1)
cell_list=find(gamma_truth(cell_group_list{this_plane}));
hold on;

for j=1:length(cell_list)
   i_cell=cell_list(j);
   trials_temp = cumsum(cell_trial_numbers(:,i_cell));
  line(trials_temp,...
             cell_gamma_mean(:,i_cell),...
            'Color','r',...
            'LineStyle','-','LineWidth',3)
end
%color_list(cell_assignments(:,i_cell),:)

%-----------------------------------%
%% Visualize the change of cell assignments 

plot(total_trial_numbers,group_assignments(:,1,j),'Color','red','LineStyle','-','LineWidth',3)

hold on;
line([15 80],[50 50],'Color','red','LineStyle','-','LineWidth',3)
text(85,50,'Number of disconnected cells')

line(n_total_trial,n_alive_cells,'Color','blue','LineStyle','-','LineWidth',3)
line([15 80],[43 43],'Color','blue','LineStyle','-','LineWidth',3)
text(85,43,'Number of connected cells')
xlim([0 max(n_total_trial)]);
ylim([0 n_cell_this_plane]);

xlabel('Number of trials');
ylabel('Cell counts');

hold off;



%%


scatter(gamma_truth,mean_gamma,'Marker','o','SizeData',25,...
    'MarkerFaceColor','b', 'MarkerEdgeColor','b', 'MarkerFaceAlpha',0.8)
x=[0 1];y=[0 1];
hold on;
line(x,y,'Color','red','LineStyle','--')
hold off;

xlim([0 1]);
ylim([0 1]);

xlabel('True synaptic success rate');
ylabel('Estimated synaptic success rate');
% saveas(1,strcat('./Figures/Toy/','Gamma_truth_vs_fits','.jpg'));


    %%
     figure(2)
    scatter(gain_truth(find(gamma_truth>0)),mean_gain(find(gamma_truth>0)),'Marker','o','SizeData',25,...
                'MarkerFaceColor','g', 'MarkerEdgeColor','g', 'MarkerFaceAlpha',0.8)
            x=[0 1];y=[0 1];
            hold on;
    line(x,y,'Color','red','LineStyle','--')
    hold off;
    
xlim([0.01 0.03]);
ylim([0.01 0.03]);
    
    xlabel('True optical gain');
ylabel('Estimated optical gain');
% saveas(2,strcat('./Figures/Toy/','Gain_truth_vs_fits','.jpg'));

    %% Number of trials v.s. dead cells and alive cells 
    n_dead_cells = zeros(final_iter,1);
    n_alive_cells = zeros(final_iter,1);
    n_trial_proc=zeros(final_iter,1);
    for i=2:final_iter
       n_dead_cells(i)=sum(dead_cells{i}); 
        n_alive_cells(i)=sum(alive_cells{i});
        n_trial_proc(i)=length(mpp_undefined{i-1})+length(mpp_connected{i-1})+length(mpp_disconnected{i-1});
    end
    n_total_trial=cumsum(n_trial_proc);
    
      %%
     figure(3)
    
     plot(n_total_trial,n_dead_cells,'Color','red','LineStyle','-','LineWidth',3)
    
            hold on;
     line([15 80],[50 50],'Color','red','LineStyle','-','LineWidth',3)
    text(85,50,'Number of disconnected cells')

     line(n_total_trial,n_alive_cells,'Color','blue','LineStyle','-','LineWidth',3)
     line([15 80],[43 43],'Color','blue','LineStyle','-','LineWidth',3)
    text(85,43,'Number of connected cells')
    xlim([0 max(n_total_trial)]);
    ylim([0 n_cell_this_plane]);
    
    xlabel('Number of trials');
    ylabel('Cell counts');
   
     hold off;
    
% saveas(3,strcat('./Figures/Toy/','Cell_counts_vs_ntrial','.jpg'));

    
    %% Plot the change over course 
    for iter = 1:final_iter
            figure(iter)
            scatter(cell_locations(cell_group_list{this_plane},2),...
                cell_locations(cell_group_list{this_plane},1),...
                'Marker','o','SizeData',1,...
                'MarkerFaceColor','r', 'MarkerEdgeColor','r', 'MarkerFaceAlpha',0.2)
            hold on;
            for i_cell = 1:length(cell_group_list{this_plane})
                cell_index =cell_group_list{this_plane}(i_cell);
                
                    if undefined_cells{iter}(i_cell)==1
                        facecolor='k';markershape='o';
                    elseif dead_cells{iter}(i_cell)==1
                        facecolor='r';markershape='x';
                    elseif alive_cells{iter}(i_cell)==1
                      facecolor='b';markershape='o';
                    elseif potentially_connected_cells{iter}(i_cell)==1
                      facecolor='g';markershape='o';
                    elseif potentially_disconnected_cells{iter}(i_cell)==1
                      facecolor='r';markershape='o';
                    end
                    
                    scatter(cell_locations(cell_index,2),...
                        cell_locations(cell_index,1),...
                        'Marker',markershape,'SizeData',200,...
                        'MarkerFaceColor',facecolor, 'MarkerEdgeColor',facecolor,...
                        'MarkerFaceAlpha',0.6)
                  
                
                hold on;
            end
            xlabel(strcat(num2str(n_total_trial(iter)), ' trials'));
            xlim([-160 160]);
ylim([-160 160]);

           %xlabel('X (um)');
            %ylabel('Y (um)');
%             axis off;
            
      
 txt=cell([5 1]);
txt{1}='Undefined';
txt{2}='Disconnected';
txt{3}='Connected';
txt{4}='Pot. connected';
txt{5}='Pot. disconnected';
xcord=[-150 -150 -150 -150 -150];
ycord=[150 130 110 90 70];
for i=1:5
    if i==1
        facecolor='k';markershape='o';
    elseif i==2
        facecolor='r';markershape='x';
    elseif i==3
        facecolor='b';markershape='o';
    elseif i==4
        facecolor='g';markershape='o';
    elseif i==5
        facecolor='r';markershape='o';
    end
scatter(xcord(i),ycord(i),...
    'Marker',markershape,'SizeData',200,...
                        'MarkerFaceColor',facecolor, 'MarkerEdgeColor',facecolor,...
                        'MarkerFaceAlpha',0.6)
                    text(xcord(i)+10,ycord(i),txt{i})
end 
    hold off;

% saveas(iter,strcat('./Figures/Toy/','Cell_map', num2str(iter),'.jpg'));

end
