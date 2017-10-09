addpath(genpath('../../../mapping-inference'));
%% Parameters:
d=15; % distance between the two cells 
%% Cell locations
cell_locations=zeros(1,3);
n_cell=1;
%% Cell parameters 
background_rate=1e-4;
v_th_known=15*ones([n_cell,1]);v_reset_known=-1e4*ones([n_cell,1]);
g_truth = 0.02*ones([n_cell,1]);
funcs.invlink = @invlink_sig;%@(resp) log(1 + exp(resp));%@(x) exp(x);
gamma_truth = [0.7];gain_truth=[0.015];
%%  Load the current template
load('../Environments/chrome-template-3ms.mat');
downsamp=1;max_time=300;power_level = 30:10:100;
num_power_level=length(power_level);
current_template=template(1:downsamp:max_time);
t_vect= 1:1:max_time;
time_max=max_time;
%% Preprocessing
% Calculate the firing probability 
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
eff_stim_threshold=stim_grid(min(find(prob_trace>0.99)));
%% Select the stimulation locations
% Load the shape template
load('../Environments/l23_template_cell.mat');
%l23_average_shape
temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;shape_template=l23_average_shape;
r1=5;r2=10;r3=15;num_per_grid=12;
num_per_grid_dense=16;
% The list of all related cells :
grid_jitters = zeros(num_per_grid,2);
for i_grid = 1:num_per_grid
   grid_jitters(i_grid,:)=[sin(2*pi*i_grid/num_per_grid) cos(2*pi*i_grid/num_per_grid)]; 
end
grid_jitters=[grid_jitters zeros(num_per_grid,1)];

% Calculate the stimulation locations 
target_locations = zeros(n_cell*(2*num_per_grid+1),3);
for i_cell=1:n_cell
    nucleus_loc=cell_locations(i_cell,:);
    grid_locs=nucleus_loc;
    grid_locs=[grid_locs;bsxfun(@plus,grid_jitters*r1,nucleus_loc)];
    grid_locs=[grid_locs;bsxfun(@plus,grid_jitters*r2,nucleus_loc)];
    target_idx=(i_cell-1)*(2*num_per_grid+1) +(1: (2*num_per_grid+1));
    target_locations(target_idx,:) = grid_locs;
end
% target_locations(:,3)= mean(cell_locations(target_cell_list.primary,3));

%plot(target_locations{this_plane}(:,2),target_locations{this_plane}(:,1),'.')


cell_params.locations =  cell_locations;
cell_params.shape_gain = ones(n_cell,1);
cell_template = struct();
cell_template.shape= shape_template;
% [pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
%     cell_template,target_locations);
[pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
    cell_template,target_locations);
%power_selected = power_level(1)*ones([size(target_locations,1) 1]);
loc_to_cell = zeros(size(pi_target,2),1);
loc_to_cell(1:(2*num_per_grid+1))=1;
loc_to_cell((2*num_per_grid+1)+ (1:(2*num_per_grid+1)))=2;
power_selected=zeros(n_cell*(2*num_per_grid+1),1);
power_sd=zeros(n_cell*(2*num_per_grid+1),1);

%% Parameters in the design stage
% Design parameters
n_spots_per_trial = 1;
% Need to run sims to check how these parameters affect the results
% Prior distribution
prior_pi0=0.8;

% Initialize tuning parameters in the VI
C_threshold = 0.01;maxit=1000;
S=200;epsilon=0.01;eta_logit=0;
background_rt=background_rate*time_max; % raw probability of firing within a trial
eta_beta=0.05;

gamma_estimates = 0.5*ones(n_cell,1);% for drawing samples (not really used)
prob_weight=0;
id_continue=1;% an indicator

% lklh_func=@calculate_likelihood_sum_bernoulli; % likelihood function is
% specificed when fitting the working model

stim_threshold = 10;
% bounds of the gamma:
gain_bound.up=0.03;gain_bound.low=0.005;


% Initialize storage
mean_gamma_current=zeros(n_cell,1);
mean_gain_current=gain_template*ones(n_cell,1);
gamma_path=zeros(n_cell,1);
gain_path=zeros(n_cell,1);
var_gamma_path=zeros(n_cell,1);

%% Simulate the prior distribution of gain:

var_gain_prior=1;

K=50;
n_replicates=1;
cell_list= 1:2;
trial_max=4000;
iter=1;
n_MC_samples=100;
epislon=1;
connected_threshold=0.4;
%% 
run('Twoneurons_naive.m');
%%
run('Twoneurons_optimal.m');
%% Plot:
figure(1)

hold on;

line((1: iter)*K*n_cell,...
    gamma_path_naive(1,2:end),...
    'Color',[1 0 0 0.5],'LineStyle','-','LineWidth',3)
% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gamma_25_path_naive(1,2:end),...
% 'Color',[1 0 0 0.2], 'LineStyle','-','LineWidth',1)
% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gamma_75_path_naive(1,2:end),...
%     'Color',[1 0 0 0.2],  'LineStyle','-','LineWidth',1)

line((1: iter)*K*n_cell,...
    gamma_path_naive(2,2:end),...
    'Color',[1 0 0 0.5],...
    'LineStyle',':','LineWidth',3)
% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gamma_25_path_naive(2,2:end),...
% 'Color',[1 0 0 0.2], 'LineStyle',':','LineWidth',1)
% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gamma_75_path_naive(2,2:end),...
%     'Color',[1 0 0 0.2],  'LineStyle',':','LineWidth',1)


% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gamma_path_optimal(1,2:end),...
%     'Color',[0 0 0 0.5],...
%     'LineStyle','-','LineWidth',3)
% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gamma_25_path_optimal(1,2:end),...
% 'Color',[0 0 0 0.2], 'LineStyle','-','LineWidth',1)
% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gamma_75_path_optimal(1,2:end),...
%     'Color',[0 0 0 0.2],  'LineStyle','-','LineWidth',1)
% 
% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gamma_path_optimal(2,2:end),...
%     'Color',[0 0 0 0.5],...
%     'LineStyle',':','LineWidth',3)
% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gamma_25_path_optimal(2,2:end),...
% 'Color',[0 0 0 0.2], 'LineStyle',':','LineWidth',1)
% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gamma_75_path_optimal(2,2:end),...
%     'Color',[0 0 0 0.2],  'LineStyle',':','LineWidth',1)


line([K*n_cell trial_max], [gamma_truth(1) gamma_truth(1)],...
    'Color',[0 1 0 0.5],...
    'LineStyle','-','LineWidth',3)
line([K*n_cell trial_max], [gamma_truth(2) gamma_truth(2)],...
    'Color',[0 1 0 0.5],...
    'LineStyle',':','LineWidth',3)

xlabel('Number of trials','FontSize',15);
ylabel('E[\gamma|data]','FontSize',25);
hold off;
%%
figure(2)

hold on;

line((1:iter)*K*n_cell,...
    gain_path_naive(1,2:end),...
    'Color',[1 0 0 0.5],...
    'LineStyle','-','LineWidth',3)

% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gain_25_path_naive(1,2:end),...
% 'Color',[1 0 0 0.2], 'LineStyle','-','LineWidth',1)
% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gain_75_path_naive(1,2:end),...
%     'Color',[1 0 0 0.2],  'LineStyle','-','LineWidth',1)

line((1: iter)*K*n_cell,...
    gain_path_naive(2,2:end),...
    'Color',[1 0 0 0.5],...
    'LineStyle',':','LineWidth',3)

% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gain_25_path_naive(2,2:end),...
% 'Color',[1 0 0 0.2], 'LineStyle',':','LineWidth',1)
% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gain_75_path_naive(2,2:end),...
%     'Color',[1 0 0 0.2],  'LineStyle',':','LineWidth',1)



% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gain_path_optimal(1,2:end),...
%     'Color',[0 0 0 0.5],...
%     'LineStyle','-','LineWidth',3)
% 
% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gain_25_path_optimal(1,2:end),...
% 'Color',[0 0 0 0.2], 'LineStyle','-','LineWidth',1)
% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gain_75_path_optimal(1,2:end),...
%     'Color',[0 0 0 0.2],  'LineStyle','-','LineWidth',1)
% 
% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gain_path_optimal(2,2:end),...
%     'Color',[0 0 0 0.5],...
%     'LineStyle',':','LineWidth',3)
% 
% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gain_25_path_optimal(2,2:end),...
% 'Color',[0 0 0 0.2], 'LineStyle',':','LineWidth',1)
% line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%     gain_75_path_optimal(2,2:end),...
%     'Color',[0 0 0 0.2],  'LineStyle',':','LineWidth',1)


line([K*n_cell trial_max], [gain_truth(1) gain_truth(1)],...
    'Color',[0 1 0 0.5],...
    'LineStyle','-','LineWidth',3)

line([K*n_cell trial_max], [gain_truth(2) gain_truth(2)],...
    'Color',[0 1 0 0.5],...
    'LineStyle',':','LineWidth',3)

xlabel('Number of trials','FontSize',15);
ylabel('E[\phi|data]','FontSize',25);
hold off;
