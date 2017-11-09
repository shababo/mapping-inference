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
eff_stim_threshold=stim_grid(min(find(prob_trace>1e-1)));
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
% n_cell_this_plane=10;
% target_cell_list(1).primary=cell_group_list{this_plane}(1:10);
target_cell_list(1).secondary=[];
related_cell_list=[target_cell_list.primary; target_cell_list.secondary]; 

v_th_known_related= v_th_known(related_cell_list);
v_reset_known_related=v_reset_known(related_cell_list);
g_related = g_truth(related_cell_list);

gamma_truth(related_cell_list(2))=0.9;
gamma_truth(related_cell_list(6))=0.9;

gamma_related = gamma_truth(related_cell_list);
% gamma_truth(3)=0.8; gamma_truth(10)=0.8;gamma_truth(1)=0.8;gamma_truth(5)=0.8;
gain_related=gain_truth(related_cell_list);

%% Select the stimulation locations
% Load the shape template
% load('./Environments/l23_template_cell.mat');
% load('./Environments/l23_cells_for_sim.mat');
%l23_average_shape
temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;shape_template=l23_average_shape;
r1=5;r2=10;r3=15;num_per_grid=12;
num_per_grid_dense=16;
stim_threshold=eff_stim_threshold/gain_template;
grid_type = 1; % 1: circles; 2: lines;

target_locations_selected=cell_locations(target_cell_list.primary,:);
target_locations_selected(:,3)=mean(target_locations_selected(:,3));

cell_params.locations =  cell_locations(related_cell_list,:);
cell_params.shape_gain = ones(length(related_cell_list),1);
cell_template = struct();
cell_template.shape= shape_template;
[pi_target_selected, inner_normalized_products] = get_weights_v2(cell_params, ...
    cell_template,target_locations_selected);

power_selected = power_level(3)*ones(length(target_cell_list.primary),1);
%num_z_grid=8;% 6 values on z-plane (1 N N N N 1)
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
K_undefined=80; % each cell appears approximately 10*2 times
K_disconnected=5; % each cell appears approximately 10*2 times
K_connected=10; % each cell appears approximately 10*2 times


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
%% Conduct random single-spot trials on the cell nuclei
mpp_undefined{iter}=[];trials_locations_undefined{iter}=[];trials_powers_undefined{iter}=[];

    cell_list= find(undefined_cells{iter});
    gamma_estimates = 0.5*ones(length(cell_list),1);% for drawing samples...
    [trials_locations, trials_powers] = random_design(target_locations_selected,power_selected,...
        inner_normalized_products,single_spot_threshold,gamma_estimates,prob_weight,...
        id_notconnected, loc_to_cell,cell_list,n_spots_per_trial,K_undefined,n_replicates);
    trials_powers=reshape(randsample(30:10:100,size(trials_powers,1)*n_spots_per_trial,true),...
        [size(trials_powers,1) n_spots_per_trial]);
    
    [cells_probabilities_undefined, stim_size_undefined] = get_prob_and_size(...
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


%% Online design:
% Initialize the path of variational families with infor from previous
% iteratin
variational_params_path.pi(:,iter+1)=var_pi_ini*ones(length(related_cell_list),1);
variational_params_path.alpha(:,iter+1)=variational_params_path.alpha(:,iter);
variational_params_path.beta(:,iter+1)=variational_params_path.beta(:,iter);
variational_params_path.alpha_gain(:,iter+1)=variational_params_path.alpha_gain(:,iter);
variational_params_path.beta_gain(:,iter+1)=variational_params_path.beta_gain(:,iter);


%------------------------------------------------------%
% Fit VI on Group A: the undefined cells
% Fit the VI on group C: potentially connected cells
% This step is different, we shoul fit each neuron seperately if possible
mean_gamma_undefined=zeros(length(related_cell_list),1);
variance_gamma_undefined=ones(length(related_cell_list),1);

cell_list= find(undefined_cells{iter});
% designs_remained=stim_size_connected(:,cell_list);

active_trials=find(sum(stim_size_undefined(:,cell_list),2)>stim_threshold);
neighbour_list= 1:length(related_cell_list);
    
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
    
    designs_remained=stim_size_undefined(active_trials,neighbour_list);
    %             active_trials=find(sum(designs_remained,2)>stim_threshold);
    %             designs_remained=designs_remained(active_trials,:);
    mpp_remained=mpp_undefined{iter}(active_trials);
    
    lklh_func=@calculate_likelihood_bernoulli;
    %             lklh_func=@lif_glm_firstevent_loglikelihood_for_VI;
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
    [mean_gain_temp, ~] = calculate_posterior_mean(...
        parameter_history.alpha_gain(:,end),parameter_history.beta_gain(:,end),gain_bound.low,gain_bound.up);
%%
%     gamma_related-mean_gamma_temp
% gain_related(neighbour_list)-mean_gain_temp
gamma_related(gamma_related>0)-mean_gamma_temp(gamma_related>0)
abs(gain_related(gamma_related>0)-mean_gain_temp(gamma_related>0))./gain_related(gamma_related>0)


%%
mean_gain_temp(gamma_related==0)
mean_gamma_temp(gamma_related==0)






