%% Load the data set for cell locations
load('../Environments/6_3_s2c3_mpp_and_stim_data.mat')
cell_locations=cell_locs;
n_cell = size(cell_locations,1);
clear mpp;
clear stim_pow;
clear target_inds;
%% Generate the cell parameters
background_rate=1e-4;
v_th_known=15*ones([n_cell,1]);
v_reset_known=-1e4*ones([n_cell,1]);
g_truth = 0.02*ones([n_cell,1]);

funcs.invlink = @invlink_sig;%@(resp) log(1 + exp(resp));%@(x) exp(x);
time_max =300;
%gamma_truth=zeros(n_cell_this_plane,1);
gamma_truth = (rand([n_cell 1])<0.2).*(0.7+0.3*rand([n_cell 1]));
% gamma_truth(3)=0.8; gamma_truth(10)=0.8;gamma_truth(1)=0.8;gamma_truth(5)=0.8;
gain_truth=0.015+rand([n_cell 1])*0.01;

if gain_type == 2 %change the gain distribution for connected cell to have a long left tail. 
    gain_lowerbound =0.008;
    gain_upperbound =0.02;
    
    n_connected=sum(gamma_truth>0);
    lognormal_temp=exp(normrnd(0,2,[n_cell 1]))/2;
    gain_truth= (gain_upperbound-gain_lowerbound)*exp(-lognormal_temp)./(1+exp(-lognormal_temp))+...
        gain_lowerbound;
%     gain_truth(gamma_truth>0) =
end


%%  Load the current template
load('../Environments/chrome-template-3ms.mat');
downsamp=1;max_time=300;power_level = 30:10:100;
num_power_level=length(power_level);
current_template=template(1:downsamp:max_time);
t_vect= 1:1:max_time;

%% Select a plane 
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


%% Simulate the prior distribution of gain:

% bounds of the gamma:
gain_bound.up=0.03;gain_bound.low=0.005;
var_alpha_initial=1;var_beta_initial=1.78;
if good_prior == 1 
    var_alpha_gain_initial=log( (gain_related+0.002*rand([length(gain_related), 1]) - gain_bound.low)./(gain_bound.up-gain_related));
elseif good_prior==2
var_alpha_gain_initial=log( (0.02 - gain_bound.low)./(gain_bound.up-0.02))*ones(length(related_cell_list),1);
 
end
var_beta_gain_initial=0.5; % uncertainty of the prior 

% The prioir info serves as the first variational distributions 
variational_params_path.alpha=var_alpha_initial*ones(length(related_cell_list),1);
variational_params_path.beta=var_beta_initial*ones(length(related_cell_list),1);
variational_params_path.alpha_gain=var_alpha_gain_initial;
variational_params_path.beta_gain=var_beta_gain_initial*ones(length(related_cell_list),1);


