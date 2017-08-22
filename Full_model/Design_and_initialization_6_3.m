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

cell_params.g=0.02;cell_params.v_th_known=15;cell_params.gain_template = 0.02;
stim_unique = (1:1000)/10;
linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};
 [prob_trace]=get_firing_probability(...
    linkfunc,current_template,stim_unique,cell_params,delay_params);
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
this_plane = 2;

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
gamma_truth = (rand([n_cell_this_plane 1])<0.1).*(0.5+0.5*rand([n_cell_this_plane 1]));
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
%
[pi_target_selected, inner_normalized_products,target_locations_selected,power_selected,...
    target_locations_all,weakly_identified_cells] = ...
    get_stim_locations(...
    cell_list,cell_locations,power_level,...
    r1,r2,r3,num_per_grid,num_per_grid_dense,shape_template,...
    stim_unique,prob_trace);
%% End of preprocessing
%-------------------------------------------------%

%% Designing experiment 
%-------------------------------------------------%
% In each batch,
%   a) design random trials on the remianing cells
%   b) design random trials on the potnetially disconnected cells 
%       (identified in the previous iteration) 

%% Parameters in the design stage

% Design parameters
n_spots_per_trial = 4;
n_replicates=2; % conduct two replicates for each trial
K_random=10; % each cell appears approximately 20 times 
K_confirm=20; % each cell appears approximately 10 times 
single_spot_threshold=5; % switch to single spot stimulation if there are fewer than 5 cells 
trial_max=2000;
disconnect_threshold = 0.15;
confirm_threshold = 0.1;

potentially_disconnected_cells=[];
confirmed_disconnected_cells=[];
% Prior distribution 
prior_pi0=0.7;

% Initialize storage 
cells_history = cell(0);
cells_history{1}=ones(n_cell_this_plane,1);
iter=1;
mpp_confirm=cell(0);
trials_locations_confirm=cell(0);
trials_powers_confirm=cell(0);

mpp_random=cell(0);
trials_locations_random=cell(0);
trials_powers_random=cell(0);

variational_params_path.pi=zeros(n_cell_this_plane,1);
variational_params_path.log_alpha=zeros(n_cell_this_plane,1);
variational_params_path.log_beta=zeros(n_cell_this_plane,1);

designs_random=[];designs_confirm=[];
outputs_random=[];outputs_confirm=[];

% Initialize the variational family
var_pi_ini=0.01;
var_log_alpha_initial=0;
var_log_beta_initial=0;


% Initialize the parameters in the VI
C_threshold = 0.01;maxit=1000;
S=200;epsilon=0.01;eta_logit=0;eta_beta=0.01;
background_rt=background_rate*time_max;


visualized = 1;
   
n_trials=0;
gamma_estimates = 0.5*ones(n_cell_this_plane,1);% for drawing samples...

id_continue=1;% an indicator 
prob_weight=0;
%% Online design: 
while ((n_trials < trial_max) & (id_continue>0))
    % while not exceeding the set threshold of total trials
    % and there are new cells being excluded 
    
    %---------------------------------
    % Random trials on the surviving cells
    % We focus more random trials on the cells that have low gammas 
    K=K_random;
    remaining_cell_list= setdiff(find(cells_history{iter}),potentially_disconnected_cells);
    
    [trials_locations,  trials_powers] = random_design(...
        target_locations_selected,power_selected,...
        inner_normalized_products,single_spot_threshold,...
        gamma_estimates,prob_weight,...
        remaining_cell_list,n_spots_per_trial,K,n_replicates);
    [cells_probabilities_random, ~] = get_prob_and_size(...
        pi_target_selected,trials_locations,trials_powers,...
        stim_unique,prob_trace);
    
    % Generate mpp given the trials
    [mpp_temp] = draw_samples(...
        trials_locations, trials_powers, pi_target_selected, background_rate,...
        v_th_known, v_reset_known, g_truth, gain_truth,gamma_truth,...
        current_template,  funcs,    delay_params,stim_threshold,time_max);
    mpp_random{iter}=mpp_temp;
    trials_locations_random{iter}=trials_locations;
    trials_powers_random{iter}=trials_powers;
    
    %-------
    % Confirm and eliminate the disconnected cells
    % i.e., the cell-killing mode
    
    if ~isempty(potentially_disconnected_cells)
        % Find cells with close to zero gammas
        K=K_confirm;
        gamma_estimates_confirm = 0.5*ones(length(potentially_disconnected_cells),1);% for drawing samples...

        remaining_cell_list= potentially_disconnected_cells;
        [trials_locations,  trials_powers] = random_design(...
            target_locations_selected,power_selected,...
            inner_normalized_products,single_spot_threshold,...
            gamma_estimates_confirm,0,...
            remaining_cell_list,n_spots_per_trial,K,n_replicates);
        
        [cells_probabilities_confirm, ~] = get_prob_and_size(...
            pi_target_selected,trials_locations,trials_powers,...
            stim_unique,prob_trace);
        
        % Conduct trials
        [mpp_temp] = draw_samples(...
            trials_locations, trials_powers, pi_target_selected, background_rate,...
            v_th_known, v_reset_known, g_truth, gain_truth,gamma_truth,...
            current_template,  funcs,    delay_params,stim_threshold,time_max);
        mpp_confirm{iter}=mpp_temp;
        trials_locations_confirm{iter}=trials_locations;
        trials_powers_confirm{iter}=trials_powers;
    else
        mpp_confirm{iter}=[];
        trials_locations_confirm{iter}=[];
        trials_powers_confirm{iter}=[];
    end

%------------------------------------------%
% Transform the data
designs_random=[designs_random; cells_probabilities_random];
n_previous_trials_random =length(outputs_random);
for i_trial = 1:size(cells_probabilities_random,1)
    outputs_random(n_previous_trials_random+i_trial,1)=length(mpp_random{iter}(i_trial).times);
end

if ~isempty(potentially_disconnected_cells)
    designs_confirm=[designs_confirm; cells_probabilities_confirm];
    n_previous_trials_confirm =length(outputs_confirm);
    for i_trial = 1:size(cells_probabilities_confirm,1)
        outputs_confirm(n_previous_trials_confirm+i_trial,1)=length(mpp_confirm{iter}(i_trial).times);
    end
end

%------------------------------------------%
% Analysis: 

% Fit VI on the remaining cells 
remaining_cell_list= find(cells_history{iter});
n_cell_remaining=length(remaining_cell_list);
variational_params=struct([]);
if iter== 1
    for i_cell = 1:n_cell_remaining
        variational_params(i_cell).pi = var_pi_ini;
        variational_params(i_cell).p_logit = log(variational_params(i_cell).pi/(1-variational_params(i_cell).pi));
        variational_params(i_cell).log_alpha = var_log_alpha_initial; %log_alpha
        variational_params(i_cell).log_beta = var_log_beta_initial;%log_alpha
    end
       prior_params.pi0= prior_pi0*ones(n_cell_remaining,1);
    prior_params.alpha0= ones(n_cell_remaining,1);
    prior_params.beta0 = ones(n_cell_remaining,1);
 
else
    for i_cell_idx = 1:n_cell_remaining
        i_cell=remaining_cell_list(i_cell_idx);
        variational_params(i_cell_idx).pi = variational_params_path.pi(i_cell,iter-1);
        variational_params(i_cell_idx).p_logit = log(variational_params(i_cell_idx).pi/(1-variational_params(i_cell_idx).pi));
        variational_params(i_cell_idx).log_alpha = variational_params_path.log_alpha(i_cell,iter-1);
        variational_params(i_cell_idx).log_beta = variational_params_path.log_beta(i_cell,iter-1);
    end 
    prior_params.pi0= [variational_params(:).pi]';
    prior_params.alpha0= exp([variational_params(:).log_alpha]');
    prior_params.beta0 = exp([variational_params(:).log_beta]');
end

% Include only the remaining cells 

designs_remained=cells_probabilities_random;
designs_remained=designs_remained(:,remaining_cell_list);
active_trials=find(sum(designs_remained,2)>1e-3);
designs_remained=designs_remained(active_trials,:);
outputs_remained=outputs_random( (n_previous_trials_random+1) : length(outputs_random));
outputs_remained=outputs_remained(active_trials,:);
%
[parameter_history,~] = fit_working_model_vi(...
    designs_remained,outputs_remained,background_rt, ...
     variational_params,prior_params,C_threshold,...
    S,epsilon,eta_logit,eta_beta,maxit);

% Record the variational parameters
variational_params_path.pi(remaining_cell_list,iter) = parameter_history.pi(:,end);
variational_params_path.log_alpha(remaining_cell_list,iter) = log(parameter_history.alpha(:,end));
variational_params_path.log_beta(remaining_cell_list,iter) = log(parameter_history.beta(:,end));

last_iter = size(parameter_history.pi,2);
mean_gamma_temp= (1-parameter_history.pi(:,last_iter)).*...
     (C_threshold+ (1-C_threshold)./(1+parameter_history.beta(:,last_iter)./parameter_history.alpha(:,last_iter)));   
%variance_gamma = (1-parameter_history.pi(:,last_iter)).* ...
%    (mean_gamma_random.^2+ parameter_history.alpha(:,last_iter).*parameter_history.beta(:,last_iter)./...
%    (parameter_history.alpha(:,last_iter)+parameter_history.beta(:,last_iter)).^2./...
%    (parameter_history.alpha(:,last_iter)+parameter_history.beta(:,last_iter)+1));
mean_gamma_random=zeros(n_cell_this_plane,1);
mean_gamma_random(remaining_cell_list,1)=mean_gamma_temp;
    


% Fit the VI on the confirmatory data sets
if ~isempty(potentially_disconnected_cells)
    remaining_cell_list= potentially_disconnected_cells;
    n_cell_remaining=length(potentially_disconnected_cells);
    variational_params=struct([]);
    if iter== 1
        for i_cell = 1:n_cell_remaining
            variational_params(i_cell).pi = var_pi_ini;
            variational_params(i_cell).p_logit = log(variational_params(i_cell).pi/(1-variational_params(i_cell).pi));
            variational_params(i_cell).log_alpha = var_log_alpha_initial; %log_alpha
            variational_params(i_cell).log_beta = var_log_beta_initial;%log_alpha
        end
    else
        for i_cell_idx = 1:n_cell_remaining
            i_cell=remaining_cell_list(i_cell_idx);
            variational_params(i_cell_idx).pi = variational_params_path.pi(i_cell,iter-1);
            variational_params(i_cell_idx).p_logit = log(variational_params(i_cell_idx).pi/(1-variational_params(i_cell_idx).pi));
            variational_params(i_cell_idx).log_alpha = variational_params_path.log_alpha(i_cell,iter-1);
            variational_params(i_cell_idx).log_beta = variational_params_path.log_beta(i_cell,iter-1);
        end
    end
    prior_params.pi0= prior_pi0*ones(n_cell_remaining,1);
    prior_params.alpha0= ones(n_cell_remaining,1);
    prior_params.beta0 = ones(n_cell_remaining,1);
    
    % Include only the remaining cells
    
designs_remained=cells_probabilities_confirm;
designs_remained=designs_remained(:,remaining_cell_list);
active_trials=find(sum(designs_remained,2)>1e-3);
designs_remained=designs_remained(active_trials,:);
outputs_remained=outputs_confirm( (n_previous_trials_confirm+1) : length(outputs_confirm));
outputs_remained=outputs_remained(active_trials,:);
   
    %
    [parameter_history,~] = fit_working_model_vi(...
        designs_remained,outputs_remained,background_rt, ...
        variational_params,prior_params,C_threshold,...
        S,epsilon,eta_logit,eta_beta,maxit);
    
    % Eliminate cells if their new estimates are still lower than the threshold
    
last_iter = size(parameter_history.pi,2);
    mean_gamma_temp= (1-parameter_history.pi(:,last_iter)).*...
        (C_threshold+ (1-C_threshold)./(1+parameter_history.beta(:,last_iter)./parameter_history.alpha(:,last_iter)));
    mean_gamma_confirm=zeros(n_cell_this_plane,1);
    mean_gamma_confirm(remaining_cell_list,1)=mean_gamma_temp;
    confirmed_disconnected_cells = find( mean_gamma_confirm<confirm_threshold );
    confirmed_disconnected_cells = intersect(confirmed_disconnected_cells,potentially_disconnected_cells);
else
    confirmed_disconnected_cells=[];
end

iter=iter+1;
cells_history{iter}=cells_history{iter-1};
cells_history{iter}(confirmed_disconnected_cells)=0;

potentially_disconnected_cells =find(mean_gamma_random<disconnect_threshold);
potentially_disconnected_cells =intersect(potentially_disconnected_cells,find(cells_history{iter}));

n_trials = size(designs_random,1)+size(designs_confirm,1);
remaining_cell_list= setdiff(find(cells_history{iter}),potentially_disconnected_cells);
    gamma_estimates=mean_gamma_random(remaining_cell_list);
    
    if isempty(potentially_disconnected_cells)
        if id_continue==2
            id_continue=0;% terminate
        else
            id_continue=2;
        end
    else
       id_continue=1; 
    end
% Plot the progress

    fprintf('Number of trials so far: %d; number of cells killed: %d\n',n_trials, sum(cells_history{iter}==0))

    if visualized == 1
        pi=variational_params_path.pi(:,iter-1);
        alpha=exp(variational_params_path.log_alpha(:,iter-1));
        beta=exp(variational_params_path.log_beta(:,iter-1));
        mean_gamma_temp= (1-pi).*...
            (C_threshold+ (1-C_threshold)./(1+beta./alpha));
        variance_gamma_temp = (1-pi).* ...
            (mean_gamma_temp.^2+ alpha.*beta./...
            (alpha+beta).^2./(alpha+beta+1));
        mean_gamma_temp(find(cells_history{iter}==0))=0.01;
        variance_gamma_temp(find(cells_history{iter}==0))=0.01;
        
        
        figure(iter)
        scatter(cell_locations(cell_group_list{this_plane},2),...
            cell_locations(cell_group_list{this_plane},1),...
            'Marker','o','SizeData',1,...
            'MarkerFaceColor','r', 'MarkerEdgeColor','r', 'MarkerFaceAlpha',1)
        hold on;
        for i_cell = 1:length(cell_group_list{this_plane})
            cell_index =cell_group_list{this_plane}(i_cell);
            
            if cells_history{iter}(i_cell)>0
            if sum(potentially_disconnected_cells==i_cell)>0
                facecolor='g';
            else
                facecolor='k';
            end
            scatter(cell_locations(cell_index,2),...
                cell_locations(cell_index,1),...
                'Marker','o','SizeData',- 200/log(variance_gamma_temp(i_cell)),...
                'MarkerFaceColor',facecolor, 'MarkerEdgeColor',facecolor,...
                'MarkerFaceAlpha',0.3+ 0.7*mean_gamma_temp(i_cell))
            else 
            % the cell has been eliminated 
            facecolor='r';
            scatter(cell_locations(cell_index,2),...
                cell_locations(cell_index,1),...
                'Marker','x','SizeData',80,...
                'MarkerFaceColor',facecolor, 'MarkerEdgeColor',facecolor)
            
            end
            hold on;
        end
        
        hold off;
        %xlabel('X (um)');
        %ylabel('Y (um)');
        axis off;
        
    end



end
%% End of the cell-killing mode
%--------------------------------------------------------%
%% Check posterior means:
% mean_gamma(gamma_truth==0)
% mean_gamma(gamma_truth>0)
% gamma_truth(gamma_truth>0)
% Visualize the posterior distributions 

%% Plot the selected cells 
for i_figure= 2:iter
    saveas(i_figure,strcat('./Figures/Preprocess/Killed_cells',num2str(i_figure),'.jpg'));
end
%% Get the numbers 
summ_matrix = zeros(size(cells_history,2),3);
% remaining cells, percentage of connected cells, percentage of disconnected cells
for i_figure= 2:iter
    summ_matrix(i_figure,1)=sum(cells_history{i_figure});
    summ_matrix(i_figure,2)=sum(gamma_truth>0 & cells_history{i_figure})/sum(gamma_truth>0);
    summ_matrix(i_figure,3)=sum(gamma_truth==0 & cells_history{i_figure})/sum(gamma_truth==0);
end

%% Second stage: learning the cell properties
%---------------------------------------------%
