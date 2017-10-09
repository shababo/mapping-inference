addpath(genpath('../../../mapping-inference'));
%% Parameters:
d=5; % distance between the two cells
%% Cell locations
cell_locations=zeros(2,3);
cell_locations(2,1)=d;
n_cell=2;
%% Cell parameters
num_sim=2;
num_seed=20;

background_rate=1e-4;
v_th_known=15*ones([n_cell,1]);v_reset_known=-1e4*ones([n_cell,1]);
g_truth = 0.02*ones([n_cell,1]);
funcs.invlink = @invlink_sig;%@(resp) log(1 + exp(resp));%@(x) exp(x);
gamma_truth = [0.7; 0.6];gain_truth=[0.005; 0.02];

%%  Load the current template
load('../Environments/chrome-template-3ms.mat');
downsamp=1;max_time=300;power_level = 30:2:100;
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
eff_stim_threshold=stim_grid(min(find(prob_trace>0.01)));
fire_stim_threshold=stim_grid(min(find(prob_trace>0.99)));

%% Select the stimulation locations
% Load the shape template
load('../Environments/l23_template_cell.mat');
%l23_average_shape
temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;shape_template=l23_average_shape;
r1=5;r2=10;num_per_grid=12;
num_per_grid_dense=16;
% The list of all related cells :
grid_jitters = zeros(num_per_grid,2);
for i_grid = 1:num_per_grid
    grid_jitters(i_grid,:)=[sin(2*pi*i_grid/num_per_grid) 0];
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

K=5;
n_replicates=1;
cell_list= 1:2;
trial_max=200;
iter=1;
n_MC_samples=100;
epislon=10;
connected_threshold=0.4;
%%

% Initialize storage
mean_gamma_current=zeros(n_cell,1);
mean_gain_current=gain_template*ones(n_cell,1);
gamma_path=zeros(n_cell,1);
gain_path=zeros(n_cell,1);
var_gamma_path=zeros(n_cell,1);
var_gain_path=zeros(n_cell,1);

gamma_25_path=zeros(n_cell,1);
gamma_75_path=zeros(n_cell,1);

gain_25_path=zeros(n_cell,1);
gain_75_path=zeros(n_cell,1);


% bounds of the gamma:
gain_bound.up=0.03;gain_bound.low=0.005;
var_alpha_initial=1;var_beta_initial=1.78;
var_alpha_gain_initial=log( (gain_truth- gain_bound.low)./(gain_bound.up-gain_truth));
var_beta_gain_initial=var_gain_prior; % uncertainty of the prior
% The prioir info serves as the first variational distributions
variational_params_path.alpha=var_alpha_initial*ones(n_cell,1);
variational_params_path.beta=var_beta_initial*ones(n_cell,1);
variational_params_path.alpha_gain=var_alpha_gain_initial;
variational_params_path.beta_gain=var_beta_gain_initial*ones(n_cell,1);

n_trials = 0;
n_events=0;
iter=1;

mpp_naive_optimal=cell(0);
trials_locations_optimal=cell(0);
trials_powers_optimal=cell(0);

%%
%---- Select stim locations and power:
% draw samples of gains from the posterior distribution
var_gain_prior_list=[0.01 0.1 0.5 1];
for i_var = 1:length(var_gain_prior_list)


var_gain_prior=var_gain_prior_list(i_var);
% bounds of the gamma:
gain_bound.up=0.03;gain_bound.low=0.005;
var_alpha_initial=1;var_beta_initial=1.78;
var_alpha_gain_initial=log( (gain_truth- gain_bound.low)./(gain_bound.up-gain_truth));
var_beta_gain_initial=var_gain_prior; % uncertainty of the prior
% The prioir info serves as the first variational distributions
variational_params_path.alpha=var_alpha_initial*ones(n_cell,1);
variational_params_path.beta=var_beta_initial*ones(n_cell,1);
variational_params_path.alpha_gain=var_alpha_gain_initial;
variational_params_path.beta_gain=var_beta_gain_initial*ones(n_cell,1);

gain_samples=zeros(n_MC_samples,n_cell);
for i_cell = 1:n_cell
    v_alpha_gain = variational_params_path.alpha_gain(i_cell,iter);
    v_beta_gain = exp(variational_params_path.beta_gain(i_cell,iter));
    temp=normrnd(v_alpha_gain,v_beta_gain,[n_MC_samples 1]);
    gain_samples(:,i_cell) = exp(temp)./(1+exp(temp))*(gain_bound.up-gain_bound.low) +gain_bound.low;
end

firing_prob=zeros(size(pi_target,2),length(power_level),size(pi_target,1));
for i_loc = 1:size(pi_target,2)
    for k=1:length(power_level)
        for i_cell = 1:n_cell
            stimulation_received=pi_target(i_cell,i_loc)*power_level(k);
            
            effective_stim= stimulation_received*gain_samples(:,i_cell);
            stim_index=max(1,round(effective_stim*stim_scale));
            prob_collapsed=sum(prob_trace_full(stim_index,:),2);
            firing_prob(i_loc,k,i_cell)=mean(prob_collapsed);
        end
    end
end
target_location_optimal=zeros(2,3);
loc_to_cell_optimal=1:n_cell;

firing_prob_difference= firing_prob(:,:,1)-...
    firing_prob(:,:,2);

% Plotting 
[x_sort,index_sort] = sort(target_locations(:,1));
% flip_index = length(power_level)+1-(1: length(power_level));
flip_index=1: length(power_level);
fp_sorted=firing_prob_difference(index_sort,flip_index);
 figure(1)
hold on;
hm=heatmap(fp_sorted');
colorbar
xlabel('X coordinate','FontSize',25);
ylabel('Power level','FontSize',25);
xticks([1 10 20 30 40 50])
xticklabels({string(round(x_sort([1 10 20 30 40 50]),1))})
yticks([1 12 24 36])   
xlim([1 50]);
ylim([1 length(power_level)]);

yticklabels({string(round(power_level([1 12 24 36]),1))})
title(strcat('var(\phi;data)=',num2str(var_gain_prior)),'FontSize',25)

scatter([ min(find(x_sort>0)) min(find(x_sort>d))],...
         [1 1],...
            'Marker','o','SizeData',50,...
            'MarkerFaceColor',[1 0 0], 'MarkerEdgeColor','w',...
            'MarkerFaceAlpha',1)

hold off;
saveas(1,strcat('./Figures/Oct03/ObjectiveVariance',num2str(i_var),'d',num2str(d),'Weak.png'));

close all;
end
