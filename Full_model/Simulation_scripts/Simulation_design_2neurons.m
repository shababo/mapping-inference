addpath(genpath('../../../mapping-inference'));
%% Specify the setting in this simulation

% for i_sim = [12 15]
%     for i_seed = 1:20

cell_parameters_type=1;
prior_info_type=1;
model_type=1;
design_type =1;
assignment_type=1;
%% Generate cellular parameters
rng(i_seed,'twister');
d=10;
background_rate=1e-4;

gamma_truth =[0.7 0];
gain_truth=[0.02 0.007];
cell_params=struct([]);

cell_locations=[0 0 0; d 0 0 ];
funcs.invlink = @invlink_sig;%@(resp) log(1 + exp(resp));%@(x) exp(x);

for i_cell = 1:2
    cell_params(i_cell).v_thresh= 15;
    cell_params(i_cell).v_reset= -1e4;
    cell_params(i_cell).location=cell_locations(i_cell,:);
    cell_params(i_cell).g=0.02;
    cell_params(i_cell).optical_gain=gain_truth(i_cell);
    cell_params(i_cell).gamma=gamma_truth(i_cell);
    cell_params(i_cell).shape_truth=1; % for simulation 
    cell_params(i_cell).shape=1; 
end
%% Preprocessing
load('../Environments/chrome-template-3ms.mat');
downsamp=1;time_max=300;power_level = 30:10:100;
current_template=template(1:downsamp:time_max);
t_vect= 1:1:time_max;
%% Specify delay distributions
delay_params.type=2; %1: normal; 2: gamma
delay_params.mean=58; delay_params.std=15;
delay_params.delayed=true; delay_params.n_grid=200;
% bounds of the gamma:
gain_bound.up=0.03;gain_bound.low=0.003;
gamma_bound.up=1;gamma_bound.low=0.01;
%% Pre-calculate the first spike intensity as a function of actual stimulation 
template_cell.optical_gain = 1; % for calculation
template_cell.g=0.02;template_cell.v_thresh=15;

max_actual_stimulation=5;
num_stim_grid=1000;
linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};

gain_template=0.02;
stim_scale=num_stim_grid/max_actual_stimulation;
stim_grid = (1:num_stim_grid)/stim_scale; % grid of actual stimulation 

% stim_unique=(1:1000)/stim_scale/gain_template;
[prob_trace_full,~] = get_first_spike_intensity(...
    linkfunc,current_template,stim_grid,template_cell,delay_params);
prob_trace=sum(prob_trace_full,2);

% Minimum actual stimulation in order for a cell to fire 
minimum_stim_threshold=stim_grid(min(find(prob_trace>0.01))); 
% The threshold of actual stimulation when a cell fires with high
% probability 
fire_stim_threshold=stim_grid(min(find(prob_trace>0.99)));
stim_threshold = minimum_stim_threshold/gain_bound.up;
%% Pre-calculate the candidate stimulation locations 
target_cell_list=struct;
target_cell_list(1).primary=[1 2];
target_cell_list(1).secondary=[];

load('../Environments/l23_template_cell.mat');
% load('../Environments/l23_cells_for_sim.mat');

temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;
shape_template=l23_average_shape;

grid_params=struct([]);
grid_params(1).radius=[5];
grid_params(1).number=[12];
grid_params(2).radius=[10];
grid_params(2).number=[16];


[pi_target, inner_normalized_products,target_locations,loc_to_cell,...
    target_locations_nuclei,pi_target_nuclei, loc_to_cell_nuclei]= ...
   get_stim_locations(target_cell_list,cell_params,grid_params,shape_template);


%% For simulation:
sim_indicator=true;
[pi_target_truth, ~,~,~,...
    ~,pi_target_nuclei_truth,~]= ...
   get_stim_locations(target_cell_list,cell_params,grid_params,shape_template,[],[],sim_indicator);

%% End of preprocessing
%-------------------------------------------------%


%% Designing experiment
%-------------------------------------------------%
%% Parameters in the design stage
n_replicates=1; % number of replicates for each trial
n_spots_per_trial = 1;trial_max=200;
K_connected=5;
% determine how many batches are used in fitting the model, default is zero 
num_trace_back = 2;
% Define tuning parameters in the VI
maxit=1000;S=50;epsilon=0.01;eta=1;eta_max=2;
background_rt=background_rate*time_max; n_MC_samples=50;

%%
related_cell_list=[target_cell_list.primary target_cell_list.secondary];
n_related_cell=length(related_cell_list);
%% Initialize the posterior/variational distributions
switch prior_info_type
    case 1 
        var_gain_prior=0;
        gain_bias=0;
        disp('Use good prior')
    case 2 %   
        var_gain_prior=0; 
        gain_bias=0;
        disp('Use uninformative prior')
    case 3 %
        var_gain_prior=0;
        gain_bias=0.005;
        disp('Biased prior')
    otherwise
end

prior_pi0=0.5;var_pi_ini=0.5;
var_alpha_initial=0;var_beta_initial=1;
perturbed_gain=max(0.0055,gain_truth(related_cell_list)+gain_bias*(2*(rand([n_related_cell, 1])-0.5) ));
var_alpha_gain_initial=log((perturbed_gain  - gain_bound.low)./(gain_bound.up-perturbed_gain));
var_beta_gain_initial=var_gain_prior; % uncertainty of the prior

prior_params=struct;
temp=num2cell(log(prior_pi0/(1-prior_pi0))*ones(n_related_cell,1));[prior_params(1:n_related_cell).p_logit]=temp{:};
temp=num2cell(var_alpha_initial*ones(n_related_cell,1));[prior_params(1:n_related_cell).alpha]=temp{:};
temp=num2cell(var_beta_initial*ones(n_related_cell,1));[prior_params(1:n_related_cell).beta]=temp{:};
temp=num2cell(var_alpha_gain_initial);[prior_params(1:n_related_cell).alpha_gain]=temp{:};
temp=num2cell(var_beta_gain_initial*ones(n_related_cell,1));[prior_params(1:n_related_cell).beta_gain]=temp{:};

variational_params_path=struct;
variational_params_path=prior_params;

quantiles_prob=[0.05 0.95];
% Initialize storage for the fitted parameters in the experiment
[parameter_path] = calculate_posterior(prior_params,gamma_bound,gain_bound,quantiles_prob,true);
%% Online design:
model_type=1; % 1 working, 2 ss 3 fuparamell
design_type_single=1; %0 random 1 optimal 

%% Online design:
iter=1;
n_trials=0;
tic
tstart=toc;


% Initialize the trial records 
mpp_history=struct([]);
connected_cells=cell(0);
connected_cells{1}=[1 1];
%%
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
     num_trials_connected=sum(connected_cells{iter})*K_connected;
         
   
    %-------
    % Conduct trials on group C, the potentially connected cells
    if num_trials_connected>0
        % Find cells with close to zero gammas
        cell_list= find(connected_cells{iter});
        [trials_locations,  trials_powers] = random_design(...
                 num_trials_connected, target_locations_nuclei,power_level,loc_to_cell_nuclei,...
                 related_cell_list,cell_list,pi_target_nuclei, [],...
            variational_params,n_MC_samples,gamma_bound,gain_bound,prob_trace_full,...
            fire_stim_threshold,stim_scale, ...
            1,K_connected,n_replicates,[],[],[],[],0,design_type_single);
       
       [mpp_temp] = draw_samples(...
            trials_locations, trials_powers, pi_target_nuclei_truth, background_rate,...
            cell_params(related_cell_list), current_template, funcs, delay_params,stim_threshold,time_max);
        [mpp_temp.batch]=deal(iter);
        if isempty(mpp_history)
            mpp_history=mpp_temp;
        else
            mpp_history(end+(1:+length(mpp_temp)))=mpp_temp;
        end
        n_trials=n_trials+length(mpp_temp);
    end
    
    %------------------------------------------%
    % Analysis:
    
    % Initialize the path of variational families with infor from previous
    % iteratin
    variational_params_path(iter+1,:)=variational_params_path(iter,:);
    parameter_path(iter+1,:)=parameter_path(iter,:);
    
   
  
    %----------------------------------------------%
    % Fit the VI on group C: potentially connected cells
    % This step is different, we shoul fit each neuron seperately if possible
    
    if sum(connected_cells{iter})>0
        indicators_all = find(ismember([mpp_history(:).batch],iter-(0:num_trace_back) ));
        mpp_all=mpp_history(indicators_all);
        trials_locations=reshape([mpp_all(:).locations],1,[])';
        trials_powers=reshape([mpp_all(:).power],1,[])';
        stim_all = get_stim_size(pi_target_nuclei,trials_locations,trials_powers);
        cell_list= find(connected_cells{iter});
       
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
    
    
    
  
    connected_cells{iter+1}=connected_cells{iter};
    iter=iter+1;
end
tend=toc;
total_computing_time= tend-tstart;

%%
gamma_related=gamma_truth(related_cell_list);
gain_related=gain_truth(related_cell_list);
