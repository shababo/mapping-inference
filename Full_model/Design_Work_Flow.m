addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Load the data set for cell locations
load('./Environments/6_3_s2c3_mpp_and_stim_data.mat')
cell_locations=cell_locs;
n_cell = size(cell_locations,1);
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

%% Consider the second plane in this analysis
this_plane =4;
% plane 4 
rng(12242,'twister');
n_cell_this_plane = length(cell_group_list{this_plane});
background_rate=1e-4;
v_th_known=15*ones([n_cell_this_plane,1]);
v_reset_known=-1e4*ones([n_cell_this_plane,1]);
g_truth = 0.02*ones([n_cell_this_plane,1]);

funcs.invlink = @invlink_sig;%@(resp) log(1 + exp(resp));%@(x) exp(x);
stim_threshold=10;
time_max =300;
%gamma_truth=zeros(n_cell_this_plane,1);
gamma_truth = (rand([n_cell_this_plane 1])<0.1).*(0.7+0.3*rand([n_cell_this_plane 1]));
% gamma_truth(3)=0.8; gamma_truth(10)=0.8;gamma_truth(1)=0.8;gamma_truth(5)=0.8;
gain_truth=0.015+rand([n_cell_this_plane 1])*0.01;


%% Select the stimulation locations
% Load the shape template
load('./Environments/l23_template_cell.mat');
load('./Environments/l23_cells_for_sim.mat');
%l23_average_shape
temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;shape_template=l23_average_shape;

cell_list=cell_group_list{this_plane};
r1=5;r2=10;r3=15;num_per_grid=12;

num_per_grid_dense=16;
stim_threshold=eff_stim_threshold/gain_template;
%
[pi_target_selected, inner_normalized_products,target_locations_selected,power_selected,...
    target_locations_all,cell_neighbours,...
    target_locations_nuclei, power_nuclei,pi_target_nuclei, loc_to_cell_nuclei] = ...
    get_stim_locations(...
    cell_list,cell_locations,power_level,...
    r1,r2,r3,num_per_grid,num_per_grid_dense,shape_template,...
    stim_unique,prob_trace,stim_threshold);
%num_z_grid=8;% 6 values on z-plane (1 N N N N 1)
%%

%% End of preprocessing
%-------------------------------------------------%

%% Designing experiment
%-------------------------------------------------%
%% Parameters in the design stage

% Design parameters
n_spots_per_trial = 4;
n_replicates=2; % conduct two replicates for each trial
K_undefined=3; % each cell appears approximately 10*2 times
K_disconnected=3; % each cell appears approximately 10*2 times
K_connected=10; % each cell appears approximately 10*2 times

single_spot_threshold=15; % switch to single spot stimulation if there are fewer than 8 cells in this group
trial_max=2000;
disconnected_threshold = 0.2;
disconnected_confirm_threshold = 0.2;


connected_threshold = 0.5;
connected_confirm_threshold = 0.5;



% Initialize the five cell groups
undefined_cells= cell(0); undefined_cells{1}=ones(n_cell_this_plane,1);%A
potentially_disconnected_cells= cell(0); potentially_disconnected_cells{1}=zeros(n_cell_this_plane,1);%B
dead_cells= cell(0); dead_cells{1}=zeros(n_cell_this_plane,1);%D
potentially_connected_cells= cell(0); potentially_connected_cells{1}=zeros(n_cell_this_plane,1);%C
alive_cells= cell(0);
alive_cells{1}=zeros(n_cell_this_plane,1);%E

%weakly_identified_cells

% Prior distribution
prior_pi0=0.8;


iter=1;
mpp_undefined=cell(0);
trials_locations_undefined=cell(0);
trials_powers_undefined=cell(0);

mpp_disconnected=cell(0);
trials_locations_disconnected=cell(0);
trials_powers_disconnected=cell(0);

mpp_connected=cell(0);
trials_locations_connected=cell(0);
trials_powers_connected=cell(0);

designs_undefined=[];designs_connected=[];designs_disconnected=[];
outputs_undefined=[];outputs_connected=[];outputs_disconnected=[];


% Initialize the variational family
var_pi_ini=0.01;
var_log_alpha_initial=0;
var_log_beta_initial=0;
var_log_alpha_gain_initial=0;
var_log_beta_gain_initial=0;

variational_params_path.pi=var_pi_ini*ones(n_cell_this_plane,1);
variational_params_path.log_alpha=var_log_alpha_initial*ones(n_cell_this_plane,1);
variational_params_path.log_beta=var_log_alpha_initial*ones(n_cell_this_plane,1);
variational_params_path.log_alpha_gain=var_log_alpha_gain_initial*ones(n_cell_this_plane,1);
variational_params_path.log_beta_gain=var_log_alpha_gain_initial*ones(n_cell_this_plane,1);

% Initialize the parameters in the VI
C_threshold = 0.01;maxit=1000;
S=200;epsilon=0.01;eta_logit=0;eta_beta=0.01;
background_rt=background_rate*time_max;


visualized = 0;

n_trials=0;
gamma_estimates = 0.5*ones(n_cell_this_plane,1);% for drawing samples...

id_continue=1;% an indicator
prob_weight=0;

lklh_func=@calculate_likelihood_sum_bernoulli;
stim_threshold = 10;
gain_bound.up=0.03;
gain_bound.low=0.01;

id_notconnected=false;
loc_to_cell = 1:size( target_locations_selected,1);

connected=true;
 %  loc_to_cell_nuclei is from get_stim_locations 
mean_gamma_current=zeros(n_cell_this_plane,1);
mean_gain_current=gain_template*ones(n_cell_this_plane,1);
change_threshold=0.05;
gamma_path=zeros(n_cell_this_plane,1);
% Online design:
while ((n_trials < trial_max) & (id_continue>0))
    % while not exceeding the set threshold of total trials
    % and there are new cells being excluded
    
    
    % Conduct random trials
    
    % On the undefined cells
    mpp_undefined{iter}=[];
    trials_locations_undefined{iter}=[];
    trials_powers_undefined{iter}=[];
    
    if sum(undefined_cells{iter})>0
        
        cell_list= find(undefined_cells{iter});
        gamma_estimates = 0.5*ones(length(cell_list),1);% for drawing samples...
        
        [trials_locations, trials_powers] = random_design(...
            target_locations_selected,power_selected,...
            inner_normalized_products,single_spot_threshold,...
            gamma_estimates,prob_weight,...
            id_notconnected, loc_to_cell,... 
            cell_list,n_spots_per_trial,K_undefined,n_replicates);
        [cells_probabilities_undefined, ~] = get_prob_and_size(...
            pi_target_selected,trials_locations,trials_powers,...
            stim_unique,prob_trace);
        
        % Generate mpp given the trials
        [mpp_temp] = draw_samples(...
            trials_locations, trials_powers, pi_target_selected, background_rate,...
            v_th_known, v_reset_known, g_truth, gain_truth,gamma_truth,...
            current_template, funcs, delay_params,stim_threshold,time_max);
        mpp_undefined{iter}=mpp_temp;
        trials_locations_undefined{iter}=trials_locations;
        trials_powers_undefined{iter}=trials_powers;
    end
    %-------
    % Conduct trials on group B, the potentially disconnected cells
    mpp_disconnected{iter}=[];
    trials_locations_disconnected{iter}=[];
    trials_powers_disconnected{iter}=[];
    if sum(potentially_disconnected_cells{iter})>0
        % Find cells with close to zero gammas
        cell_list= find(potentially_disconnected_cells{iter});
        gamma_estimates_confirm = 0.5*ones(length(cell_list),1);% for drawing samples...
        [trials_locations,  trials_powers] = random_design(...
            target_locations_selected,power_selected,...
            inner_normalized_products,single_spot_threshold,...
            gamma_estimates_confirm,0,...
             id_notconnected, loc_to_cell,... 
            cell_list,n_spots_per_trial,K_disconnected,n_replicates);
        [cells_probabilities_disconnected, ~] = get_prob_and_size(...
            pi_target_selected,trials_locations,trials_powers,...
            stim_unique,prob_trace);
        
        % Conduct trials
        [mpp_temp] = draw_samples(...
            trials_locations, trials_powers, pi_target_selected, background_rate,...
            v_th_known, v_reset_known, g_truth, gain_truth,gamma_truth,...
            current_template,  funcs,    delay_params,stim_threshold,time_max);
        mpp_disconnected{iter}=mpp_temp;
        trials_locations_disconnected{iter}=trials_locations;
        trials_powers_disconnected{iter}=trials_powers;
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
        [trials_locations,  trials_powers] = random_design(...
            target_locations_nuclei,power_nuclei,...
            inner_normalized_products,single_spot_threshold,...
            gamma_estimates_confirm,0,...
            connected,  loc_to_cell_nuclei,... 
            cell_list,1,K_connected,n_replicates);
        %[cells_probabilities_connected, ~] = get_prob_and_size(...
        %    pi_target_nuclei,trials_locations,trials_powers,...
        %    stim_unique,prob_trace);
        [~, stim_size_connected] = get_prob_and_size(...
            pi_target_nuclei,trials_locations,trials_powers,...
            stim_unique,prob_trace);
        
        % Conduct trials
        [mpp_temp] = draw_samples(...
            trials_locations, trials_powers, pi_target_nuclei, background_rate,...
            v_th_known, v_reset_known, g_truth, gain_truth,gamma_truth,...
            current_template,  funcs,    delay_params,stim_threshold,time_max);
        mpp_connected{iter}=mpp_temp;
        trials_locations_connected{iter}=trials_locations;
        trials_powers_connected{iter}=trials_powers;
    end
    
    
    %------------------------------------------%
    % Transform the data
    % no need to record the probabilities all the time..
    
    %cells_probabilities_undefined;
    if sum(undefined_cells{iter})>0
        for i_trial = 1:size(cells_probabilities_undefined,1)
            outputs_undefined(i_trial,1)=length(mpp_undefined{iter}(i_trial).times);
        end
        n_trials=n_trials+i_trial;
    end
    if  sum(potentially_disconnected_cells{iter})>0
        %cells_probabilities_disconnected;
        for i_trial = 1:size(cells_probabilities_disconnected,1)
            outputs_disconnected(i_trial,1)=length(mpp_disconnected{iter}(i_trial).times);
        end
        n_trials=n_trials+i_trial;
    end
    if  sum(potentially_connected_cells{iter})>0
        %cells_probabilities_disconnected;
%         for i_trial = 1:size(stim_size_connected,1)
%             outputs_connected(i_trial,1)=length(mpp_connected{iter}(i_trial).times);
%         end
        n_trials=n_trials+size(stim_size_connected,1);
    end
    
    %------------------------------------------%
    % Analysis:
    
    variational_params_path.log_alpha_gain(:,iter+1)=0;
    variational_params_path.log_alpha_gain(:,iter+1)=0;
    variational_params_path.pi(:,iter+1)=var_pi_ini*ones(n_cell_this_plane,1);
    variational_params_path.log_alpha(:,iter+1)=variational_params_path.log_alpha(:,iter);
    variational_params_path.log_beta(:,iter+1)=variational_params_path.log_beta(:,iter);
    variational_params_path.log_alpha_gain(:,iter+1)=variational_params_path.log_alpha_gain(:,iter);
    variational_params_path.log_beta_gain(:,iter+1)=variational_params_path.log_beta_gain(:,iter);
    
    
    %------------------------------------------------------%
    % Fit VI on Group A: the undefined cells
    mean_gamma_undefined=zeros(n_cell_this_plane,1);
    
    if sum(undefined_cells{iter})>0
        cell_list= find(undefined_cells{iter});
        % Update variational and prior distribution
        variational_params=struct([]);
        for i_cell_idx = 1:length(cell_list)
            i_cell=cell_list(i_cell_idx);
            variational_params(i_cell_idx).pi = variational_params_path.pi(i_cell,iter);
            variational_params(i_cell_idx).p_logit = log(variational_params(i_cell_idx).pi/(1-variational_params(i_cell_idx).pi));
            variational_params(i_cell_idx).log_alpha = variational_params_path.log_alpha(i_cell,iter);
            variational_params(i_cell_idx).log_beta = variational_params_path.log_beta(i_cell,iter);
        end
        prior_params.pi0= [variational_params(:).pi]';
        prior_params.alpha0= exp([variational_params(:).log_alpha]');
        prior_params.beta0 = exp([variational_params(:).log_beta]');
        
        
        
        designs_remained=cells_probabilities_undefined(:,cell_list);
        active_trials=find(sum(designs_remained,2)>1e-3);
        designs_remained=designs_remained(active_trials,:);
        outputs_remained=outputs_undefined(active_trials,:);
        
        % find neighbours that are not in cell_list:
        neighbour_list=find(sum(cell_neighbours(cell_list,:),1)>0)';
        neighbour_list=setdiff(neighbour_list,cell_list);
        designs_neighbours=cells_probabilities_undefined(active_trials,neighbour_list);
        gamma_neighbours=mean_gamma_current(neighbour_list);
        
        lklh_func=@calculate_likelihood_bernoulli;
        % calculate_likelihood_bernoulli for multiple events 
        [parameter_history,~] = fit_working_model_vi(...
            designs_remained,outputs_remained,background_rt, ...
            variational_params,prior_params,C_threshold,...
            designs_neighbours,gamma_neighbours,...
            S,epsilon,eta_logit,eta_beta,maxit,lklh_func);
        
        % Record the variational parameters
        variational_params_path.pi(cell_list,iter+1) = parameter_history.pi(:,end);
        variational_params_path.log_alpha(cell_list,iter+1) = log(parameter_history.alpha(:,end));
        variational_params_path.log_beta(cell_list,iter+1) = log(parameter_history.beta(:,end));
        
        mean_gamma_temp= (1-parameter_history.pi(:,end)).*...
            (C_threshold+ (1-C_threshold)./(1+parameter_history.beta(:,end)./parameter_history.alpha(:,end)));
        mean_gamma_undefined=zeros(n_cell_this_plane,1);
        mean_gamma_undefined(cell_list,1)=mean_gamma_temp;
        mean_gamma_current(cell_list)=mean_gamma_temp;
        gamma_path(cell_list,iter+1)=mean_gamma_temp;
    end
    %-------------------------------------------------------------%
    
    %----------------------------------------------------------------%
    % Fit the VI on Group B: potentially disconnected cells
    mean_gamma_disconnected=ones(n_cell_this_plane,1);
    if sum(potentially_disconnected_cells{iter})>0
        cell_list= find(potentially_disconnected_cells{iter});
        variational_params=struct([]);
        for i_cell_idx = 1:length(cell_list)
            i_cell=cell_list(i_cell_idx);
            variational_params(i_cell_idx).pi = variational_params_path.pi(i_cell,iter);
            variational_params(i_cell_idx).p_logit = log(variational_params(i_cell_idx).pi/(1-variational_params(i_cell_idx).pi));
            variational_params(i_cell_idx).log_alpha = variational_params_path.log_alpha(i_cell,iter);
            variational_params(i_cell_idx).log_beta = variational_params_path.log_beta(i_cell,iter);
        end
        
        prior_params.pi0= [variational_params(:).pi]';
        prior_params.alpha0= exp([variational_params(:).log_alpha]');
        prior_params.beta0 = exp([variational_params(:).log_beta]');
        % Include only the remaining cells
        
        designs_remained=cells_probabilities_disconnected(:,cell_list);
        active_trials=find(sum(designs_remained,2)>1e-3);
        designs_remained=designs_remained(active_trials,:);
        outputs_remained=outputs_disconnected(active_trials,:);
        
         % find neighbours that are not in cell_list:
        neighbour_list=find(sum(cell_neighbours(cell_list,:),1)>0)';
        neighbour_list=setdiff(neighbour_list,cell_list);
        designs_neighbours=cells_probabilities_disconnected(active_trials,neighbour_list);
        gamma_neighbours=mean_gamma_current(neighbour_list);
       
        lklh_func=@calculate_likelihood_bernoulli;
        [parameter_history,~] = fit_working_model_vi(...
            designs_remained,outputs_remained,background_rt, ...
            variational_params,prior_params,C_threshold,...
            designs_neighbours,gamma_neighbours,...
            S,epsilon,eta_logit,eta_beta,maxit,lklh_func);
        
        % Record the variational parameters
        variational_params_path.pi(cell_list,iter+1) = parameter_history.pi(:,end);
        variational_params_path.log_alpha(cell_list,iter+1) = log(parameter_history.alpha(:,end));
        variational_params_path.log_beta(cell_list,iter+1) = log(parameter_history.beta(:,end));
        
        
        % obtain estimates
        mean_gamma_temp= (1-parameter_history.pi(:,end)).*...
            (C_threshold+ (1-C_threshold)./(1+parameter_history.beta(:,end)./parameter_history.alpha(:,end)));
        mean_gamma_disconnected=ones(n_cell_this_plane,1);
        mean_gamma_disconnected(cell_list,1)=mean_gamma_temp;
        mean_gamma_current(cell_list)=mean_gamma_temp;
        gamma_path(cell_list,iter+1)=mean_gamma_temp;
    end
    %---------------------------------------------%
    
    %----------------------------------------------%
    % Fit the VI on group C: potentially connected cells
    % This step is different, we shoul fit each neuron seperately if possible
    mean_gamma_connected=zeros(n_cell_this_plane,1);
    variance_gamma_connected=ones(n_cell_this_plane,1);
    
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
            
            neighbour_list=find(sum(cell_neighbours(cell_list(cluster_of_cells{i_cluster}),:),1)>0)';
            variational_params=struct([]);
            for i_cell_idx = 1:length(neighbour_list)
                i_cell=neighbour_list(i_cell_idx);
                variational_params(i_cell_idx).pi = variational_params_path.pi(i_cell,iter);
                variational_params(i_cell_idx).p_logit = log(variational_params(i_cell_idx).pi/(1-variational_params(i_cell_idx).pi));
                variational_params(i_cell_idx).log_alpha = variational_params_path.log_alpha(i_cell,iter);
                variational_params(i_cell_idx).log_beta = variational_params_path.log_beta(i_cell,iter);
                variational_params(i_cell_idx).log_alpha_gain = variational_params_path.log_alpha_gain(i_cell,iter);
                variational_params(i_cell_idx).log_beta_gain = variational_params_path.log_alpha_gain(i_cell,iter);
            end
            
            prior_params.pi0= [variational_params(:).pi]';
            prior_params.alpha0= exp([variational_params(:).log_alpha]');
            prior_params.beta0 = exp([variational_params(:).log_beta]');
            prior_params.alpha0_gain=  exp([variational_params(:).log_alpha_gain]');
            prior_params.beta0_gain =exp([variational_params(:).log_beta_gain]');
            
            designs_remained=stim_size_connected(:,neighbour_list);
            active_trials=find(sum(designs_remained,2)>stim_threshold);
            designs_remained=designs_remained(active_trials,:);
            mpp_remained=mpp_connected{iter}(active_trials);
            
%             
%             % find neighbours that are not in cell_list:
%             neighbour_list=find(sum(cell_neighbours(cell_list(cluster_of_cells{i_cluster}),:),1)>0)';
%             neighbour_list=setdiff(neighbour_list,cell_list(cluster_of_cells{i_cluster}));
%             designs_neighbours=stim_size_connected(active_trials,neighbour_list);
%             gamma_neighbours=mean_gamma_current(neighbour_list);
%               gain_neighbours=mean_gain_current(neighbour_list);
%       
        designs_neighbours=[];        gamma_neighbours=[];         gain_neighbours=[];
            [parameter_history] = fit_full_model_vi(...
                designs_remained, mpp_remained, background_rate, ...
                prob_trace_full,    stim_grid,...
                stim_scale,eff_stim_threshold,gain_bound,...
                variational_params,prior_params,C_threshold,stim_threshold,...
                designs_neighbours,gamma_neighbours,gain_neighbours,...
                S,epsilon,eta_logit,eta_beta,maxit);
            
                
            %      lklh_func=@calculate_likelihood_bernoulli;
            %     [parameter_history,~] = fit_working_model_vi(...
            %             designs_remained,outputs_remained,background_rt, ...
            %             variational_params,prior_params,C_threshold,...
            %             S,epsilon,eta_logit,eta_beta,maxit,lklh_func);
            %
            
            %cell_list(cluster_of_cells{i_cluster})
            variational_params_path.pi(neighbour_list,iter+1) = parameter_history.pi(:,end);
            variational_params_path.log_alpha(neighbour_list,iter+1) = log(parameter_history.alpha(:,end));
            variational_params_path.log_beta(neighbour_list,iter+1) = log(parameter_history.beta(:,end));
            variational_params_path.log_alpha_gain(neighbour_list,iter+1) = log(parameter_history.alpha_gain(:,end));
            variational_params_path.log_beta_gain(neighbour_list,iter+1) = log(parameter_history.beta_gain(:,end));
            
            % obtain estimates
            mean_gamma_temp= (1-parameter_history.pi(:,end)).*...
                (C_threshold+ (1-C_threshold)./(1+parameter_history.beta(:,end)./parameter_history.alpha(:,end)));
           
           
            mean_gain_temp= (gain_bound.low+ (gain_bound.up-gain_bound.low)...
                ./(1+parameter_history.beta_gain(:,end)./parameter_history.alpha_gain(:,end)));
           
             % approximated variance:
%             sumab=parameter_history.alpha(:,end)+parameter_history.beta(:,end);
%             var_gamma_temp=(parameter_history.alpha(:,end).*parameter_history.beta(:,end))./...
%                 (sumab.*sumab.*(sumab+1));
%             variance_gamma_connected(neighbour_list,1)=var_gamma_temp;
%            
            % mean_gamma_temp
            %mean_gain_temp
            %
            % Needs to take gain into account (as gain and gamma are
            % unidentifiable at no responses
            
           mean_gamma_connected(neighbour_list,1)=mean_gamma_temp;
           mean_gamma_current(neighbour_list)=mean_gamma_temp;
           mean_gain_current(neighbour_list)=mean_gain_temp;
           gamma_path(neighbour_list,iter+1)=mean_gamma_temp;
           
            %length([mpp_remained.times])
        end
    end
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
        change_gamma =abs(gamma_path(:,iter+1)-gamma_path(:,iter));
        connected_to_alive = intersect(find(change_gamma<change_threshold),...
            connected_to_alive);
        
        % Eliminate the weakly identifiable pairs if they are both assign to a
        % group:
        %moved_cells = [connected_to_dead; connected_to_alive]';
        %cells_and_neighbours=find(sum(cell_neighbours(moved_cells,:),1)>0);
        %neighbours_not_included=intersect(find(potentially_connected_cells{iter}), setdiff(cells_and_neighbours,moved_cells));
        %blacklist=find(sum(cell_neighbours(neighbours_not_included,:),1)>0);
        %connected_to_dead=setdiff(connected_to_dead ,blacklist);
        %connected_to_alive=setdiff(connected_to_alive,blacklist);
        
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
        % Plot the progress
        
        fprintf('Number of trials so far: %d; number of cells killed: %d\n',n_trials, sum(dead_cells{iter}+alive_cells{iter}))
       
        
        
    end
    %%
      [find(gamma_truth>0) find(alive_cells{iter})]
  
%% ------------------------------------------ 
    % Plotting 
    
    % Column I: preprocessing
    
    % I.a a 3D plot of cells
    x=cell_locations(:,1);
    y=cell_locations(:,2);
    z=cell_locations(:,3);
    
    figure(101)
    scatter3(-y,x,z,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor','b',...
        'MarkerFaceAlpha',0.6)
    view(-30,10)
    xlim([-150 150]);
    ylim([-150 150]);
    zlim([0 150]);
    saveas(101,strcat('./Figures/work_flow/','FigureIa','.jpg'));

    %%
    % I.b the 3D plot with a plane 
    x=cell_locations(:,1);
    y=cell_locations(:,2);
    z=cell_locations(:,3);
    
    figure(102)
    scatter3(-y,x,z,...
        'MarkerEdgeColor','b',...
        'MarkerFaceColor','b',...
        'MarkerFaceAlpha',0.2)
    view(-30,10)
    xlim([-150 150]);
    ylim([-150 150]);
    zlim([45 65]);
    hold on;
    
    scatter3(-y(cell_group_list{this_plane}),x(cell_group_list{this_plane}),z(cell_group_list{this_plane}),...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','k',...
        'MarkerFaceAlpha',1)
    
    [x_coord, y_coord] = meshgrid(-150:10:150);
    z_coord= mean(cell_locations(cell_group_list{this_plane},3))*ones(size(x_coord,1),size(x_coord,2));
    surf(-y_coord,x_coord,z_coord,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','none');
    hold off;
     saveas(102,strcat('./Figures/work_flow/','FigureIb','.jpg'));

    %%
    % I.c projection of cells onto that plane 
    figure(103)
    scatter(-cell_locations(cell_group_list{this_plane},2),...
        cell_locations(cell_group_list{this_plane},1),...
        'Marker','o','SizeData',30,...
        'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.8)
    hold on;
     xlim([-150 150]);
    ylim([-150 150]);
    axis off;
    hold off;
    saveas(103,strcat('./Figures/work_flow/','FigureIc','.jpg'));

    
    %% ---------------------------------------
    % Column II: heatmaps of trials
    iter_II=4;
    
    %%
    % II.a trials for undeifned cells 
    figure(201)
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
%      target_locations_selected
 saveas(201,strcat('./Figures/work_flow/','FigureIIa','.jpg'));

    
    %%
    % II.b trials for disconnected cells
    figure(202)
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
    
    facecolor='r';
    for i_trial = 1:size(trials_locations_disconnected{iter_II},1)
        scatter(target_locations_selected(trials_locations_disconnected{iter_II}(i_trial,:),2),...
            target_locations_selected(trials_locations_disconnected{iter_II}(i_trial,:),1),...
            'Marker','d','SizeData',60,...
            'MarkerFaceColor',facecolor, 'MarkerEdgeColor',facecolor,...
            'MarkerFaceAlpha',0.2)
    end
    %      target_locations_selected
    
     saveas(202,strcat('./Figures/work_flow/','FigureIIb','.jpg'));


    %%
    % II.c trials for connected cells (3D)
 figure(203)
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
     saveas(203,strcat('./Figures/work_flow/','FigureIIc','.jpg'));

    %% ---------------------------------------
    % Column III: observed data 
    %%
    % III.a Response and Probability 
    figure(3011)
    [cells_probabilities_undefined, ~] = get_prob_and_size(...
        pi_target_selected,trials_locations_undefined{iter_II},trials_powers_undefined{iter_II},...
        stim_unique,prob_trace);
    
    for i_trial = 1:size(cells_probabilities_undefined,1)
        outputs_undefined(i_trial,1)=1*isempty(mpp_undefined{iter_II}(i_trial).times);
    end
    
        
    facecolor='k';
    for i_trial = 1:size(trials_locations_undefined{iter_II},1)
        scatter(-1,...
            i_trial,...
            'Marker','s','SizeData',20,...
            'MarkerFaceColor',facecolor, 'MarkerEdgeColor','w',...
            'MarkerFaceAlpha',max(outputs_undefined(i_trial),0.01))
        hold on;
       for i_cell = 1:size(cells_probabilities_undefined,2)
           if cells_probabilities_undefined(i_trial,i_cell)>0.01
       scatter(i_cell,...
            i_trial,...
            'Marker','s','SizeData',20,...
            'MarkerFaceColor',facecolor, 'MarkerEdgeColor','w',...
            'MarkerFaceAlpha',min(max(cells_probabilities_undefined(i_trial,i_cell),0.01),1))
           end
       end
    end
    %      target_locations_selected
      xlim([-3 i_cell]);
    ylim([0 i_trial]);
  hold off;
   saveas(3011,strcat('./Figures/work_flow/','FigureIIIa1','.jpg'));

  % Fitted values
  cell_list=find(undefined_cells{iter_II});
  alpha_current=exp(variational_params_path.log_alpha(cell_list,iter_II));
  beta_current=exp(variational_params_path.log_beta(cell_list,iter_II));
  figure(3012)
 xgrid=0:0.02:1; 
    facecolor='k';
       for i_cell = 1:length(cell_list)
           
       line(xgrid,betapdf(xgrid,alpha_current(i_cell),beta_current(i_cell)));
       end
    
    %      target_locations_selected
      xlim([0 1]);
    ylim([0 5]);
  hold off;
       saveas(3012,strcat('./Figures/work_flow/','FigureIIIa2','.jpg'));

    %%
    % III.b Process and Stim size 
    
 figure(3021)
      [~, stim_size_connected] = get_prob_and_size(...
            pi_target_nuclei,trials_locations_connected{iter_II},trials_powers_connected{iter_II},...
            stim_unique,prob_trace);
           
        transparencies= stim_size_connected/max(max(stim_size_connected));
        cell_list=find(sum(transparencies)>4);
    facecolor='k';
    for i_trial = 1:size(trials_locations_undefined{iter_II},1)
        if ~isempty(mpp_connected{iter_II}(i_trial).times)
           for i_event = 1:length(mpp_connected{iter_II}(i_trial).times)
               x_sp=[mpp_connected{iter_II}(i_trial).times(i_event) mpp_connected{iter_II}(i_trial).times(i_event)];
               x_sp=x_sp/(time_max/10);
               y_sp=[i_trial+0.2 i_trial+0.8];
               line(x_sp,y_sp,'LineWidth',5)
               
           end
        end
          
        
       for i_cell = 1:length(cell_list)
           i_cell_idx = cell_list(i_cell);
           %if transparencies(i_trial,i_cell)>0.2
       scatter(i_cell+10,...
            i_trial,...
            'Marker','s','SizeData',20,...
            'MarkerFaceColor',facecolor, 'MarkerEdgeColor','w',...
            'MarkerFaceAlpha', transparencies(i_trial,i_cell_idx))
       
        hold on;
           %end
       end
    end
      xlim([-1 i_cell+10]);
    ylim([0 i_trial]);
  hold off;
    gamma_truth(find(potentially_connected_cells{iter_II}));
       saveas(3021,strcat('./Figures/work_flow/','FigureIIIb1','.jpg'));

    %% Draw the estimated gamma and gains 
      figure(3022)
  % Fitted values
  cell_list=find(potentially_connected_cells{iter_II});
  alpha_current=exp(variational_params_path.log_alpha(cell_list,iter_II));
  beta_current=exp(variational_params_path.log_beta(cell_list,iter_II));
  
 xgrid=0:0.02:1; 
  figure(4)
    facecolor='k';
       for i_cell = 1:length(cell_list)
           
       line(xgrid,betapdf(xgrid,alpha_current(i_cell),beta_current(i_cell)));
       end
    
    %      target_locations_selected
      xlim([0 1]);
    ylim([0 5]);
  hold off;
       saveas(3022,strcat('./Figures/work_flow/','FigureIIIb2','.jpg'));

  %% Gains 
  alpha_gain_current=exp(variational_params_path.log_alpha(cell_list,iter_II));
  beta_gain_current=exp(variational_params_path.log_beta(cell_list,iter_II));
  
  figure(3023)
  xgrid=0:0.02:1;
  scaled_xgrid=(gain_bound.low+ (gain_bound.up-gain_bound.low)*xgrid);
  figure(5)
  facecolor='k';
  for i_cell = 1:length(cell_list)
      
      line(scaled_xgrid,betapdf(xgrid,alpha_current(i_cell),beta_current(i_cell)));
      
  end
  
  %      target_locations_selected
  xlim([0.01 0.03]);
  ylim([0 5]);
  hold off;
     saveas(3023,strcat('./Figures/work_flow/','FigureIIIb3','.jpg'));

           
  
    %% ---------------------------------------
    % Column IV: assignments
    
    
    % IV.a undefined and potentially disconnected cells  
    cell_list=find(undefined_cells{iter_II});
    alpha_current=exp(variational_params_path.log_alpha(cell_list,iter_II));
    beta_current=exp(variational_params_path.log_beta(cell_list,iter_II));
    
    gamma_mean = alpha_current./(alpha_current+beta_current);
    xgrid=0:0.02:1;
    figure(401)
    for i_cell = 1:length(cell_list)
        if gamma_mean(i_cell)>connected_threshold
            color_codes = [0 1 0 0.4]; %green
        elseif gamma_mean(i_cell)<disconnected_threshold
            color_codes = [1 0 0 0.4]; %red
        else
            color_codes = [0 0 0 0.2]; %black 
        end
        line([gamma_mean(i_cell) gamma_mean(i_cell)],[4 5],'LineStyle','-','LineWidth',1,'Color',color_codes)
        line(xgrid,betapdf(xgrid,alpha_current(i_cell),beta_current(i_cell)),'Color',color_codes );
        hold on;
    end
    line([connected_threshold connected_threshold],[3 5],'LineStyle','-','LineWidth',1,'Color',[0 1 0])
    line([disconnected_threshold disconnected_threshold],[3 5],'LineStyle','-','LineWidth',1,'Color',[1 0 0])
    
    %      target_locations_selected
    xlim([0 1]);
    ylim([0 5]);
    hold off;
       saveas(401,strcat('./Figures/work_flow/','FigureIVa','.jpg'));

  %%
    % IV.b Potentially connected cells 
      % Fitted values
  cell_list=find(potentially_connected_cells{iter_II});
  alpha_current=exp(variational_params_path.log_alpha(cell_list,iter_II));
  beta_current=exp(variational_params_path.log_beta(cell_list,iter_II));
  gamma_mean = alpha_current./(alpha_current+beta_current);
  
   %cell_list_old=find(potentially_connected_cells{iter_II-1});
  alpha_current_old=exp(variational_params_path.log_alpha(cell_list,iter_II-1));
  beta_current_old=exp(variational_params_path.log_beta(cell_list,iter_II-1));
  gamma_mean_old = alpha_current_old./(alpha_current_old+beta_current_old);
  
  gamma_change = abs(gamma_mean-gamma_mean_old);
 change_threshold=0.06;
  xgrid=0:0.02:1; 
  figure(402)
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
       line(xgrid,betapdf(xgrid,alpha_current(i_cell),beta_current(i_cell)),'Color',color_codes );
       line(xgrid,betapdf(xgrid,alpha_current_old(i_cell),beta_current_old(i_cell)),...
           'LineStyle','--','Color',color_codes );
       
       hold on;
       end
    line([connected_threshold connected_threshold],[3 5],'LineStyle','-','LineWidth',1,'Color',[0 1 0])
    line([disconnected_threshold disconnected_threshold],[3 5],'LineStyle','-','LineWidth',1,'Color',[1 0 0])
    
    %      target_locations_selected
      xlim([0 1]);
    ylim([0 5]);
  hold off;
       saveas(402,strcat('./Figures/work_flow/','FigureIVb','.jpg'));

    %% --------------------------------------
    % Column V: new trials 
    
    iter_V=iter_II+1;
    
    %%
    % V.a trials for undeifned cells 
    figure(501)
           for i_cell = 1:length(cell_group_list{this_plane})
                cell_index =cell_group_list{this_plane}(i_cell);
                
                    if undefined_cells{iter_V}(i_cell)==1
                        facecolor='k';markershape='o';markeralpha=0.4;
                    elseif dead_cells{iter_V}(i_cell)==1
                        facecolor='r';markershape='x';markeralpha=0.2;
                    elseif alive_cells{iter_V}(i_cell)==1
                      facecolor='b';markershape='o';markeralpha=0.2;
                    elseif potentially_connected_cells{iter_V}(i_cell)==1
                      facecolor='g';markershape='o';markeralpha=0.2;
                    elseif potentially_disconnected_cells{iter_V}(i_cell)==1
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
    for i_trial = 1:size(trials_locations_undefined{iter_V},1)
           scatter(target_locations_selected(trials_locations_undefined{iter_V}(i_trial,:),2),...
                        target_locations_selected(trials_locations_undefined{iter_V}(i_trial,:),1),...
                        'Marker','d','SizeData',60,...
                        'MarkerFaceColor',facecolor, 'MarkerEdgeColor',facecolor,...
                        'MarkerFaceAlpha',0.2)
    end
%      target_locations_selected
   saveas(501,strcat('./Figures/work_flow/','FigureVa','.jpg'));

    
    %%
    % II.b trials for disconnected cells
    figure(502)
    for i_cell = 1:length(cell_group_list{this_plane})
        cell_index =cell_group_list{this_plane}(i_cell);
        
        if undefined_cells{iter_V}(i_cell)==1
            facecolor='k';markershape='o';markeralpha=0.2;
        elseif dead_cells{iter_V}(i_cell)==1
            facecolor='r';markershape='x';markeralpha=0.2;
        elseif alive_cells{iter_V}(i_cell)==1
            facecolor='b';markershape='o';markeralpha=0.2;
        elseif potentially_connected_cells{iter_V}(i_cell)==1
            facecolor='g';markershape='o';markeralpha=0.2;
        elseif potentially_disconnected_cells{iter_V}(i_cell)==1
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
    
    facecolor='r';
    for i_trial = 1:size(trials_locations_disconnected{iter_V},1)
        scatter(target_locations_selected(trials_locations_disconnected{iter_V}(i_trial,:),2),...
            target_locations_selected(trials_locations_disconnected{iter_V}(i_trial,:),1),...
            'Marker','d','SizeData',60,...
            'MarkerFaceColor',facecolor, 'MarkerEdgeColor',facecolor,...
            'MarkerFaceAlpha',0.2)
    end
    %      target_locations_selected
    
       saveas(502,strcat('./Figures/work_flow/','FigureVb','.jpg'));


    %%
    % II.c trials for connected cells (3D)
 figure(503)
    for i_cell = 1:length(cell_group_list{this_plane})
        cell_index =cell_group_list{this_plane}(i_cell);
        
        if undefined_cells{iter_V}(i_cell)==1
            facecolor='k';markershape='o';markeralpha=0.2;
        elseif dead_cells{iter_V}(i_cell)==1
            facecolor='r';markershape='x';markeralpha=0.2;
        elseif alive_cells{iter_V}(i_cell)==1
            facecolor='b';markershape='o';markeralpha=0.2;
        elseif potentially_connected_cells{iter_V}(i_cell)==1
            facecolor='g';markershape='o';markeralpha=0.5;
        elseif potentially_disconnected_cells{iter_V}(i_cell)==1
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
    for i_trial = 1:size(trials_locations_connected{iter_V},1)
           scatter3(-target_locations_nuclei(trials_locations_connected{iter_V}(i_trial,:),2),...
                        target_locations_nuclei(trials_locations_connected{iter_V}(i_trial,:),1),...
                        target_locations_nuclei(trials_locations_connected{iter_V}(i_trial,:),3),...
                        'Marker','d','SizeData',50,...
                        'MarkerFaceColor',facecolor, 'MarkerEdgeColor',facecolor,...
                        'MarkerFaceAlpha',0.2)
    end
    %xlabel(strcat(num2str(n_total_trial(iter)), ' trials'));
    view(-30,40)
    xlim([-150 150]);
    ylim([-150 150]);
    zlim([45 65]); 
       saveas(503,strcat('./Figures/work_flow/','FigureVc','.jpg'));
