addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Load the data set 
load('./Environments/6_3_s2c3_mpp_and_stim_data.mat')
target_locations =target_locs;
single_trial_limit=max(find(isnan(target_inds(:,2))));
%load('./Environments/6_3_cell_locs_reduced_and_fixed.mat')

%% Select the multi-spot data sets 
cell_locations = cell_locs;
n_cell=size(cell_locations,1);
trials_locations= target_inds( (single_trial_limit+1):end,:);
trials_powers = stim_pow((single_trial_limit+1):end);
n_trial = size(trials_locations,1);
min_time=0;max_time=300;
mpp_multi=struct([]);
for i_trial = 1:n_trial
    mpp_multi(i_trial).locations = trials_locations(i_trial,:);
    mpp_multi(i_trial).power = trials_powers(i_trial,:);
    mpp_temp=mpp(i_trial+single_trial_limit);
    if mpp_temp.num_events >0
        range_idx = mpp_temp.times<max_time & mpp_temp.times>min_time ;
        mpp_multi(i_trial).num_events = sum(range_idx);
        mpp_multi(i_trial).times = mpp_temp.times(range_idx);
        mpp_multi(i_trial).amp =mpp_temp.amp(range_idx);
    end
end
%% First step: filter the extreme cases (definitely-connected or definitely-disconnected)

% Connect the cell locations with target locations (simply use the closest
% one?)
target_to_cell = zeros(size(target_locations,1),1);
for i_loc = 1:size(target_locations,1)
    dist_temp=sum( (cell_locations(:,1:2)- target_locations(i_loc,1:2)).^2,2);
    [~, i_cell] = min(dist_temp);
    target_to_cell(i_loc)=i_cell;
end

% One issue: the closest cell in 3D might not be the closest cell in 2D 
% Calculate the number of events and number of trials for each cell:
initial_responses = zeros(n_cell,2); % first column: events; second column: trials; 
for  i_trial = 1:n_trial
    cells_stimulated = target_to_cell(mpp_multi(i_trial).locations);
    initial_responses(cells_stimulated,2)= initial_responses(cells_stimulated,2)+1;
    if ~isempty(mpp_multi(i_trial).num_events)
    initial_responses(cells_stimulated,1)=initial_responses(cells_stimulated,1)+mpp_multi(i_trial).num_events;
    end
end

% basic summary: 
sum(initial_responses(:,2)>9) % 454 cells being stimulated no less than 10 times 

%% For definitely connected cells 
gamma_threshold_connected = 0.8;
alpha=0.1;
num_trials_threshold = 9;
def_connected = zeros(n_cell,1);
for i_cell = 1:n_cell
    if initial_responses(i_cell, 2)>num_trials_threshold
        def_connected(i_cell)= initial_responses(i_cell,1)>...
        binoinv(alpha, initial_responses(i_cell, 2),gamma_threshold_connected);
    end
end
% 
 binoinv(alpha, 10,gamma_threshold_connected) % check the threshold 
 binoinv(alpha, 20,gamma_threshold_connected) % check the threshold 
 
 sum(def_connected) %50 cells are probably connected 
%% Consider the trials that do not involve these definitely connected cells 

% Calculate the number of events and number of trials for each cell (again):
connected_cell_list = 1:n_cell;connected_cell_list=connected_cell_list(def_connected==1);
initial_responses = zeros(n_cell,2); % first column: events; second column: trials; 
for  i_trial = 1:n_trial
    cells_stimulated = target_to_cell(mpp_multi(i_trial).locations);
    temp = cells_stimulated==connected_cell_list;
    if sum(sum(temp))==0
    initial_responses(cells_stimulated,2)= initial_responses(cells_stimulated,2)+1;
    if ~isempty(mpp_multi(i_trial).num_events)
    initial_responses(cells_stimulated,1)=initial_responses(cells_stimulated,1)+mpp_multi(i_trial).num_events;
    end
    end
end

% basic summary: 
unique(initial_responses(:,2))
sum(initial_responses(:,2)>5) % 403 cells left 
%% For definitely disconnected cells 
gamma_threshold_disconnected = 0.3;
alpha=0.1;
num_trials_threshold = 5;
def_disconnected = zeros(n_cell,1);
for i_cell = 1:n_cell
    if initial_responses(i_cell, 2)>num_trials_threshold
        def_disconnected(i_cell)= initial_responses(i_cell,1)<...
        binoinv(alpha, initial_responses(i_cell, 2),gamma_threshold_disconnected);
    end
end
% 
 binoinv(alpha, 8,gamma_threshold_disconnected) % check the threshold 
 binoinv(alpha, 20,gamma_threshold_disconnected) % check the threshold 
 
 sum(def_disconnected) %only 134 cells are probably disconnected.. 

 %% To summarize:
 selected_cells = def_connected | ~def_disconnected;
 n_selected = sum(selected_cells);
 selected_cell_locations = cell_locations(selected_cells,:);
 
%% Use the new initialization strategy 
% Load the shape template 
load('./Environments/l23_template_cell.mat');
%l23_average_shape
temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;
shape_template=l23_average_shape;

%------------------------------%
% Calculate the size of stimuli
cell_params.locations =  selected_cell_locations;
cell_params.shape_gain = ones(n_selected,1);
cell_template = struct();
cell_template.shape= shape_template;
[pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
    cell_template,target_locations);
stimuli_size=zeros(n_trial,n_selected);
responses=zeros(n_trial,1);
for l = 1:n_trial
    for m = 1:size(mpp_multi(l).locations,2)
        if isnan(mpp_multi(l).locations(m))
        else
            stimuli_size(l,:) = stimuli_size(l,:)+...
                (pi_target(:,mpp_multi(l).locations(m)).*mpp_multi(l).power)';
        end
    end
    if isempty(mpp_multi(l).num_events)
    else
       responses(l)=mpp_multi(l).num_events;
    end
end

%% Run a preliminary regression
% rescale the stimuli size by the maximum element on each column 
scaled_stimuli=zeros(n_trial,n_selected);

for i_cell = 1:n_selected
    scaled_stimuli(:,i_cell)=stimuli_size(:,i_cell)/max(stimuli_size(:,i_cell));
end

gamma_initial = inv(scaled_stimuli'*scaled_stimuli)*scaled_stimuli'*responses;
gamma_initial = max(0, gamma_initial);
gamma_initial = min(1, gamma_initial);

gamma_initial_crude = (scaled_stimuli'*responses)./sum(scaled_stimuli,1)';
gamma_initial_crude = max(0, gamma_initial_crude);
gamma_initial_crude = min(1, gamma_initial_crude);
%%
bcounts=15;
figure(1) 
histogram(gamma_initial,bcounts)
xlabel('Initial gamma');
ylabel('Frequency')
saveas(1,strcat('./Figures/Gibbs/','Hist_Gamma_initial_regression','.jpg'));
 
figure(2) 
histogram(gamma_initial_crude,bcounts)
xlabel('Initial gamma');
ylabel('Frequency')
saveas(2,strcat('./Figures/Gibbs/','Hist_Gamma_initial_crude','.jpg'));
 
%% Draw maps of the two initialization schemes 
figure(3)

% 
for i = 1:n_selected
    if gamma_initial(i)>0.05
    temp2 = scatter(selected_cell_locations(i,2),-selected_cell_locations(i,1),...
    200*gamma_initial(i),'MarkerEdgeColor','r','MarkerFaceColor','r',...
    'MarkerFaceAlpha',0.3);

hold on;
    end
    
    
    if gamma_initial_crude(i)>0.05
    temp2 = scatter(selected_cell_locations(i,2),-selected_cell_locations(i,1),...
    200*gamma_initial_crude(i),'MarkerEdgeColor','g','MarkerFaceColor','g',...
    'MarkerFaceAlpha',0.3);

hold on;
    end
    
end

temp5 = scatter([-150 -90 -30],[-175 -175 -175],...
    200*[0.8 0.5 0.2]);
set(temp5,'MarkerEdgeColor','r','MarkerFaceColor','r');
hold on;
txt4 = '\gamma = 0.8';
txt5 = '\gamma = 0.5';
txt6 = '\gamma = 0.2';
text(-140,-175,txt4)
text(-80,-175,txt5)
text(-20,-175,txt6)
hold off;

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])
xlabel('X (um)');
ylabel('Y (um)');

axis off;

saveas(3,strcat('./Figures/Gibbs/','Regression_vs_initial','.jpg'));

%%
%------------------------------%
% Calculate the size of stimuli
trials_locations= target_inds;
trials_powers = stim_pow;

for i_trial = 1:n_trial_all
    mpp(i_trial).locations = trials_locations(i_trial,:);
    mpp(i_trial).power = trials_powers(i_trial,:);
    mpp_temp=mpp(i_trial);
    if mpp_temp.num_events >0
        range_idx = mpp_temp.times<max_time & mpp_temp.times>min_time ;
        mpp(i_trial).num_events = sum(range_idx);
        mpp(i_trial).times = mpp_temp.times(range_idx);
        mpp(i_trial).amp =mpp_temp.amp(range_idx);
    end
end

cell_params.locations =  selected_cell_locations;
cell_params.shape_gain = ones(n_selected,1);
cell_template = struct();
cell_template.shape= shape_template;
[pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
    cell_template,target_locations);
n_trial_all = length(mpp);
stimuli_size=zeros(n_trial_all ,n_selected);
responses=zeros(n_trial_all,1);
for l = 1:n_trial_all
    for m = 1:size(mpp(l).locations,2)
        if isnan(mpp(l).locations(m))
        else
            stimuli_size(l,:) = stimuli_size(l,:)+...
                (pi_target(:,mpp(l).locations(m)).*mpp(l).power)';
        end
    end
    if isempty(mpp(l).num_events)
    else
       responses(l)=mpp(l).num_events;
    end
end

%% Run a preliminary regression
% rescale the stimuli size by the maximum element on each column 
scaled_stimuli=zeros(n_trial_all,n_selected);

for i_cell = 1:n_selected
    scaled_stimuli(:,i_cell)=stimuli_size(:,i_cell)/max(stimuli_size(:,i_cell));
end

gamma_initial = inv(scaled_stimuli'*scaled_stimuli)*scaled_stimuli'*responses;
gamma_initial = max(0, gamma_initial);
gamma_initial = min(1, gamma_initial);

gamma_initial_crude = (scaled_stimuli'*responses)./sum(scaled_stimuli,1)';
gamma_initial_crude = max(0, gamma_initial_crude);
gamma_initial_crude = min(1, gamma_initial_crude);

%%
bcounts=15;
figure(1) 
histogram(gamma_initial,bcounts)
xlabel('Initial gamma');
ylabel('Frequency')
saveas(1,strcat('./Figures/Gibbs/','Hist_Gamma_initial_regression_all','.jpg'));
 
figure(2) 
histogram(gamma_initial_crude,bcounts)
xlabel('Initial gamma');
ylabel('Frequency')
saveas(2,strcat('./Figures/Gibbs/','Hist_Gamma_initial_crude_all','.jpg'));
 
%% Draw maps of the two initialization schemes 
figure(3)

% 
for i = 1:n_selected
    if gamma_initial(i)>0.05
    temp2 = scatter(selected_cell_locations(i,2),-selected_cell_locations(i,1),...
    200*gamma_initial(i),'MarkerEdgeColor','r','MarkerFaceColor','r',...
    'MarkerFaceAlpha',0.3);

hold on;
    end
    
    
    if gamma_initial_crude(i)>0.05
    temp2 = scatter(selected_cell_locations(i,2),-selected_cell_locations(i,1),...
    200*gamma_initial_crude(i),'MarkerEdgeColor','g','MarkerFaceColor','g',...
    'MarkerFaceAlpha',0.3);

hold on;
    end
    
end

temp5 = scatter([-150 -90 -30],[-175 -175 -175],...
    200*[0.8 0.5 0.2]);
set(temp5,'MarkerEdgeColor','r','MarkerFaceColor','r');
hold on;
txt4 = '\gamma = 0.8';
txt5 = '\gamma = 0.5';
txt6 = '\gamma = 0.2';
text(-140,-175,txt4)
text(-80,-175,txt5)
text(-20,-175,txt6)
hold off;

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])
xlabel('X (um)');
ylabel('Y (um)');

axis off;

saveas(3,strcat('./Figures/Gibbs/','Regression_vs_initial_all','.jpg'));

