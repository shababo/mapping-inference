%% Load the data set for cell locations
load('../Environments/6_3_s2c3_mpp_and_stim_data.mat')
cell_locations=cell_locs;
n_cell = size(cell_locations,1);
clear mpp stim_pow target_inds;

%% Simulation some more neurons:
switch spatial_density %
    case 1 % default, original cells
        n_extra=0;
    case 2 % simulate 50% more neurons
        n_extra=ceil(n_cell/2);
    case 3 % simulate twice as many neurons
        n_extra=ceil(n_cell);    
end


extra_locations = 2*pi*rand(n_extra,3); % only need the first two columns 
for i= 1:n_extra
    temp_cell=randsample(1:n_cell,1);
    extra_locations(i,:)=cell_locations(temp_cell,:)+ ...
        distance_neighbour*[sin(extra_locations(i,1))*sin(extra_locations(i,2)),cos(extra_locations(i,1))*sin(extra_locations(i,2)),cos(extra_locations(i,2))];
end

% for i=1:3
%     extra_locations(:,i)=extra_locations(:,i)*range(cell_locations(:,i))+min(cell_locations(:,i));
% end

cell_locations=[cell_locations;extra_locations];

n_cell = size(cell_locations,1);
%% Generate the cell parameters
background_rate=1e-4;

switch sparsity
    case 1 % very sparse
        connected_proportion=0.1;
    case 2 % normal
        connected_proportion=0.2;
    case 3 % dense
        connected_proportion=0.4;
end


funcs.invlink = @invlink_sig;%@(resp) log(1 + exp(resp));%@(x) exp(x);
switch cell_parameters_type
    case 1 % normal
        gamma_truth = (rand([n_cell 1])<connected_proportion).*(0.5+0.3*rand([n_cell 1]));
        gain_truth=0.02+(rand([n_cell 1])-0.5)*0.01;
        disp('normal gains and gammas')
    case 2 % normal gamma, extreme gains
        gain_bound_truth.up=0.03;gain_bound_truth.low=0.005;
        
        gamma_truth = (rand([n_cell 1])<connected_proportion).*(0.5+0.3*rand([n_cell 1]));
        n_connected=sum(gamma_truth>0);
        lognormal_temp=exp(normrnd(0,sqrt(1),[n_cell 1]));
        gain_truth= (gain_bound_truth.up-gain_bound_truth.low)*exp(-lognormal_temp)./(1+exp(-lognormal_temp))+...
            gain_bound_truth.low;
        gain_truth=max(gain_bound_truth.low+0.0005,gain_truth);
        disp('log-normal gains, normal gammas')
    case 3 % weak gamma, normal gains
        weak_gamma_proportion=0.5; % some gammas are small
        gamma_truth = (rand([n_cell 1])<connected_proportion).*(0.5+0.3*rand([n_cell 1]));
        n_connected=sum(gamma_truth>0);
        weak_index=randsample(find(gamma_truth>0), floor(weak_gamma_proportion*n_connected));
        gamma_truth(weak_index)=(0.15+0.1*rand([length(weak_index) 1]));
        gain_truth=0.02+(rand([n_cell 1])-0.5)*0.01;
        disp('normal gains, many weak gammas')
        
end


%%
cell_params=struct([]);
for i_cell = 1:n_cell
    cell_params(i_cell).v_thresh= 15;
    cell_params(i_cell).v_reset= -1e4;
    cell_params(i_cell).location=cell_locations(i_cell,:);
    cell_params(i_cell).g=0.02;
    cell_params(i_cell).optical_gain=gain_truth(i_cell);
    cell_params(i_cell).gamma=gamma_truth(i_cell);
    cell_params(i_cell).shape_truth=1; % for simulation
    cell_params(i_cell).shape=1;
    
end

%% Put the cells into ten non-overlapping groups by their z-coordinates

n_planes = 10;
z_thresholds = quantile(cell_locations(:,3), (1:(n_planes))*0.1);
% alternatively, the thresholds can be determined by absolute depths

for i_cell = 1:size(cell_locations,1)
    cell_params(i_cell).group= sum(cell_locations(i_cell,3)>z_thresholds)+1;
end

%% Select one plane as the primary plane
%
n_cell_this_plane=sum([cell_params(:).group]==this_plane);
target_cell_list=struct([]);
target_cell_list(1).primary=find([cell_params(:).group]==this_plane);

% include also the cells in neighbouring planes to account for their inputs
% neighbour_plane = [this_plane-1 this_plane+1];
% target_cell_list(1).secondary=find(ismember([cell_params(:).group],neighbour_plane));

target_cell_list(1).secondary=[];
%% For simulation
% related_cell_list=[target_cell_list.primary; target_cell_list.secondary];
% v_th_known_related= v_th_known(related_cell_list);
% v_reset_known_related=v_reset_known(related_cell_list);
% g_related = g_truth(related_cell_list);
% gamma_related = gamma_truth(related_cell_list);
% gain_related=gain_truth(related_cell_list);
