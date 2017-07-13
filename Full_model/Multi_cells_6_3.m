addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Load real data
load('./Environments/6_3_s2c3_mpp_and_stim_data.mat')
Z_dense =target_locs;powers_trials = stim_pow;
locations_trials=target_inds;single_trial_limit=max(find(isnan(target_inds(:,2))));
locations_trials = target_inds(1:single_trial_limit,:);

powers_trials = stim_pow(1:single_trial_limit);
load('./Environments/6_3_cell_locs_reduced_and_fixed.mat')
local_locations = cell_locs;n_trial = size(locations_trials,1);n_cell_local=size(cell_locs,1);
%%
load('./Environments/chrome-template-3ms.mat');
downsamp=1;max_time=300;
power_level = unique(powers_trials);
num_power_level=length(power_level);
current_template=template(1:downsamp:max_time);
t_vect= 1:1:max_time;
v_reset_known=-4e3;
min_time=50;
for i_trial = 1:n_trial
    if mpp(i_trial).num_events >0
        range_idx = mpp(i_trial).times<max_time & mpp(i_trial).times>min_time ;
        mpp(i_trial).num_events = sum(range_idx);
        mpp(i_trial).times = mpp(i_trial).times(range_idx);
        mpp(i_trial).amp =mpp(i_trial).amp(range_idx);
    end
end
% Estimate the background rate

mpp_copy=mpp;
mpp=mpp(1:n_trial);
% Why are there a lot of events around 500?
trial_counts = 0;
event_counts = 0;
for i_trial = 1:n_trial
    %if stim_pow(i_trial)==50
    trial_counts = trial_counts+1;
    event_counts = event_counts +sum(mpp_copy(i_trial).times>20 & mpp_copy(i_trial).times<80);
    %end
end
background_rate = event_counts/trial_counts/60;

%% Other constants:
n_grid_voltage =1000;
voltage_threshold=-50;
dt=1;
n_grid_time = max_time;
t_factor=1;
sd_range=0.1;
n_delay_grid = 200;
stimulus_threshold=1e-4;
n_stimuli_grid=10;
gap_stimuli=0.05;
V_threshold=-50;
first_only=true;
funcs.invlink = @invlink_sig;%@(resp) log(1 + exp(resp));%@(x) exp(x);

load('./Environments/l23_template_cell.mat');
%l23_average_shape
temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;
cell_params.locations = local_locations;
cell_params.shape_gain = ones(n_cell_local,1);
shape_template = struct();
shape_template.shape= l23_average_shape;
[pi_dense_local, inner_normalized_products] = get_weights_v2(cell_params, shape_template,Z_dense);


stimuli_size_local=zeros(n_trial,n_cell_local);
for l = 1:n_trial
    for m = 1:size(locations_trials,2)
        if isnan(locations_trials(l,m))
        else
            stimuli_size_local(l,:) = stimuli_size_local(l,:)+(pi_dense_local(:,locations_trials(l,m)).*powers_trials(l))';
        end
    end
end
n_trial = size(stimuli_size_local,1);


% Count the number of events in trials at each locations (so that we can
% select the trials 
unique_locs_all =unique(locations_trials(:,1));
event_prop = zeros(length(unique_locs_all),1);
for i=1:length(unique_locs_all)
   relevant_index =locations_trials(:,1)==unique_locs_all(i); 
   event_prop(i) = length([mpp(relevant_index).times])/sum(relevant_index);
end
%% Draw the cell map
%this_cell=133;
figure(11)
%temp2 = scatter(local_locations(this_cell,2)+151,-local_locations(this_cell,1)-151,...
%    100,'MarkerEdgeColor','g','MarkerFaceColor','b',...
%    'MarkerFaceAlpha',0.5);
%hold on;
colors= ['w' 'g' 'b' 'k'];
prop_thresholds = [0 0.2 0.3];
colors_ind = zeros(length(unique_locs_all),1);
for i=1:length(unique_locs_all)
   colors_ind(i) = sum(event_prop(i)>prop_thresholds)+1;
end

cell_labels = cellstr(num2str( [1:size(local_locations,1)]' ));
temp1 = scatter(local_locations(:,2)+151,-local_locations(:,1)-151,...
    40);
hold on;
set(temp1,'MarkerEdgeColor','k','MarkerFaceColor','white');
alpha(temp1,1);
text(local_locations(:,2)+151,-local_locations(:,1)-151,cell_labels,'col','b');
hold on;

for i=1:length(unique_locs_all)
temp3 = scatter(target_locs(unique_locs_all(i),2)+151,-target_locs(unique_locs_all(i),1)-151,...
    40,'MarkerFaceColor',colors(colors_ind(i)));
hold on;
end

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])
xlim([-20,313]);
ylim([-313,20]);
xlabel('X (um)');
ylabel('Y (um)');
% saveas(11,strcat('./Full_cell_map.png'));

%% Pick a cell and identify its neighbours
this_cell=134;
dist= zeros(size(local_locations,1),1);
for i= 1:size(local_locations,1)
    dist(i) = sqrt(sum( (local_locations(this_cell,:) -local_locations(i,:)).^2));
end
selected_cells = find(dist<25);
selected_cells
n_selected_cells = length(selected_cells);
%%
% selected_cells = [79 85 83 90 91];
%% Pick trials that are close to the selected cells
max_radius = 25; % pick the maximum radius to consider
num_locs = size(target_locs,1);

diffs = target_locs - local_locations(this_cell,:);
relevant_locs = sqrt(sum(diffs.^2,2))<max_radius;
temp = 1:num_locs;
relevant_locs_index = temp(relevant_locs);

stim_counts = zeros(length(relevant_locs_index),1);
relevant_trials_index = [];
for i_trial = 1:n_trial
    if sum(relevant_locs_index == locations_trials(i_trial,1) )>0
        relevant_trials_index = [relevant_trials_index i_trial];
        stim_counts(find(relevant_locs_index == locations_trials(i_trial,1)))=...
            stim_counts(find(relevant_locs_index == locations_trials(i_trial,1)))+1;
    end
end

mpp_this_cell=mpp(relevant_trials_index);
locations_trials_this_cell = locations_trials(relevant_trials_index,:);
powers_trials_this_cell = powers_trials(relevant_trials_index);
stim_this_cell = stimuli_size_local(relevant_trials_index,selected_cells);
length(relevant_trials_index)
%% Draw the responses
unique_locs= unique(locations_trials_this_cell(:,1));

figure(10)
colors=['r' 'g' 'b' 'k' 'y'];
temp2 = scatter(local_locations(selected_cells ,2)+151,-local_locations(selected_cells ,1)-151,...
    500,'MarkerEdgeColor','k','MarkerFaceColor','white',...
    'MarkerFaceAlpha',0.5);
hold on;
xlim([min(local_locations(selected_cells ,2)+151)-20,max(local_locations(selected_cells ,2)+151)+20]);
ylim([-max(local_locations(selected_cells ,1)+151)-20,-min(local_locations(selected_cells ,1)+151)+20]);
xlabel('X (um)');
ylabel('Y (um)');

% saveas(10,strcat('./Two_cells_system.png'));
%%

figure(11)
colors=['r' 'g' 'b' 'k' 'y'];
temp2 = scatter(local_locations(selected_cells ,2)+151,-local_locations(selected_cells ,1)-151,...
    500,'MarkerEdgeColor','k','MarkerFaceColor','white',...
    'MarkerFaceAlpha',0.5);
hold on;
for i = 1:length(unique_locs)
    temp1 = scatter(target_locs(unique_locs(i) ,2)+151,-target_locs(unique_locs(i) ,1)-151,...
        150,'MarkerEdgeColor','k','MarkerFaceColor',colors(i),...
        'MarkerFaceAlpha',1);
    hold on;
end
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])
xlim([min(local_locations(selected_cells ,2)+151)-10,max(local_locations(selected_cells ,2)+151)+10]);
ylim([-max(local_locations(selected_cells ,1)+151)-10,-min(local_locations(selected_cells ,1)+151)+10]);
xlabel('X (um)');
ylabel('Y (um)');
% saveas(1,strcat('./Two_cells_system.png'));

% saveas(11,strcat('./Two_cells_system_stims.png'));


%%

%
for i = 1:length(unique_locs)
    figure(i)
    mpp_this_loc = mpp_this_cell(locations_trials_this_cell(:,1)==unique_locs(i));
    powers_trials_this_loc = powers_trials_this_cell(locations_trials_this_cell(:,1)==unique_locs(i));
    
    plot(sort([mpp_this_loc(powers_trials_this_loc==100).times]),1:length([mpp_this_loc(powers_trials_this_loc==100).times]),...
        '.','col',colors(i),'MarkerSize',20)
    hold on;
    plot(sort([mpp_this_loc(powers_trials_this_loc==75).times]),1:length([mpp_this_loc(powers_trials_this_loc==75).times]),...
        '.','col',colors(i),'MarkerSize',15)
    hold on;
    plot(sort([mpp_this_loc(powers_trials_this_loc==50).times]),1:length([mpp_this_loc(powers_trials_this_loc==50).times]),...
        '.','col',colors(i),'MarkerSize',10)
    
    hold off;
    ylim([0,20]);
    xlim([0,300]);
    xlabel('Time (1/20 ms)');
    ylabel('Indices (irrelevant)');
    
%     saveas(i,strcat('./Trials', num2str(i),'.png'));
end
%end

%% Prepare the small system for analysis
% Obtain the reduced trials:

% mpp_this_cell=mpp(relevant_trials_index);
% locations_trials_this_cell = locations_trials(relevant_trials_index,:);
% powers_trials_this_cell = powers_trials(relevant_trials_index);
% stim_this_cell = stimuli_size_local(relevant_trials_index,this_cell);
% unique_locs= unique(locations_trials_this_cell(:,1));

%%
n_trial_this_cell = length(mpp_this_cell);
responses = zeros(n_trial_this_cell, length(current_template));
stims = zeros(n_selected_cells,n_trial_this_cell, length(current_template));

for i_trial = 1:n_trial_this_cell
    for i_cell = 1:n_selected_cells
        k=stim_this_cell(i_trial,i_cell);
        stim = current_template*k;
        stims(i_cell, i_trial,:) = stim;
    end
    spikes_delay =round(mpp_this_cell(i_trial).times);
    responses(i_trial,spikes_delay)=1;
end
%% Delay paramters in inferece
% delay_params.mean=60;
% delay_params.std=20;


delay_params.type=2; %1: normal; 2: gamma

delay_params.mean=100;
delay_params.std=30;
delay_params.delayed=true;
delay_params.n_grid=200;

v_reset_known= -4e3*ones(n_selected_cells,1);
v_th_known=15*ones(n_selected_cells,1);

loss_type=2;
gain_lower_bound=0.002;
gain_upper_bound=0.5;
in_params.g =  0.2*ones(n_selected_cells,1);
single_connection=true;

%% Grab some sense of the gain parameters 
too_early = round(median([mpp.times]))-delay_params.mean;
cell_params.gain = max(max(stim_this_cell));
gain_sequence = 0.002+[1:100]/5e2; %from 0.03 to 0.2
n_stimuli_grid_temp=length(gain_sequence);

cell_params.V_th=v_th_known(1);
cell_params.V_reset=v_reset_known(i_cell);
cell_params.g =0.2;
cell_params.gain_sd=0.01;
     
     
[Stimuli_grid, Intensity_grid]=Intensity_v8(...
    gain_sequence,current_template,... % data from exp
    cell_params,... % estimated parameters
    funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
    n_stimuli_grid_temp,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
    V_threshold,stimulus_threshold,first_only);

for i_grid = 1:n_stimuli_grid_temp
    if sum(Intensity_grid{i_grid}(1:too_early)) > 0.5
        gain_crude= gain_sequence(i_grid);
        break;
    end
end
%% Search the connectivity one by one  
responses_reg=responses;
stims_reg=stims;
% Only assigh values in the first column 
gain_grid =[0 0.02 0.04]+gain_crude;
initial_values = zeros(length(gain_grid),2*n_selected_cells);
initial_values(:,1)=gain_grid;
[stats_conv,pred_prob] = fit_lifglm_v6(responses_reg, stims_reg,in_params,...
    background_rate,v_reset_known,v_th_known,first_only,loss_type,...
    gain_lower_bound,gain_upper_bound,initial_values,delay_params,single_connection);
% Output:
stats_conv
sum(sum(sum(pred_prob)))
%% Draw the predicted intensity 
 unique_locs= unique(locations_trials_this_cell(:,1));
n_grid = length(current_template);
figure(10)
colors=['r' 'g' 'b' 'k'];
temp2 = scatter(local_locations(selected_cells ,2)+151,-local_locations(selected_cells ,1)-151,...
    500,'MarkerEdgeColor','k','MarkerFaceColor','white',...
    'MarkerFaceAlpha',0.5);
hold on;

temp2 = scatter(local_locations(selected_cells ,2)+151,-local_locations(selected_cells ,1)-151,...
    1200*stats_conv((1:n_selected_cells)*2)+0.1,'MarkerEdgeColor','k','MarkerFaceColor','k',...
    'MarkerFaceAlpha',0.3);
hold on;

for i = 1:length(unique_locs)
    temp1 = scatter(target_locs(unique_locs(i) ,2)+151,-target_locs(unique_locs(i) ,1)-151,...
        50,'MarkerEdgeColor','k','MarkerFaceColor',colors(i),...
        'MarkerFaceAlpha',1);
    hold on;
end
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])
xlim([min(local_locations(selected_cells ,2)+151)-3,max(local_locations(selected_cells ,2)+151)+3]);
ylim([-max(local_locations(selected_cells ,1)+151)-3,-min(local_locations(selected_cells ,1)+151)+3]);
xlabel('X (um)');
ylabel('Y (um)');
%  saveas(10,strcat('./Two_cells_system_stims_fits.png'));
for i = 1:length(unique_locs)
    figure(i)
    mpp_this_loc = mpp_this_cell(locations_trials_this_cell(:,1)==unique_locs(i));
    powers_trials_this_loc = powers_trials_this_cell(locations_trials_this_cell(:,1)==unique_locs(i));
    
    plot(sort([mpp_this_loc(powers_trials_this_loc==100).times]),1:length([mpp_this_loc(powers_trials_this_loc==100).times]),...
        '.','col',colors(i),'MarkerSize',20)
    hold on;
    plot(sort([mpp_this_loc(powers_trials_this_loc==75).times]),1:length([mpp_this_loc(powers_trials_this_loc==75).times]),...
        '.','col',colors(i),'MarkerSize',15)
    hold on;
    plot(sort([mpp_this_loc(powers_trials_this_loc==50).times]),1:length([mpp_this_loc(powers_trials_this_loc==50).times]),...
        '.','col',colors(i),'MarkerSize',10)
    hold on;
    
    i_trial = min(find(locations_trials_this_cell(:,1) == unique_locs(i)&powers_trials_this_cell==75 ) );
    prob_this_loc = pred_prob(2,i_trial,:)+background_rate;
    prob_this_loc = reshape(prob_this_loc,[n_grid 1]);
    plot(1:n_grid,prob_this_loc*100,'col',colors(i),'LineStyle','--')
    hold on;
    
     i_trial = min(find(locations_trials_this_cell(:,1) == unique_locs(i)&powers_trials_this_cell==100 ) );
    prob_this_loc = pred_prob(2,i_trial,:)+background_rate;
    prob_this_loc = reshape(prob_this_loc,[n_grid 1]);
    plot(1:n_grid,prob_this_loc*100,'col',colors(i))
    
    hold off;
    ylim([0,20]);
    xlim([0,300]);
    xlabel('Time (1/20 ms)');
    ylabel('Indices (irrelevant)');
    
%     saveas(i,strcat('./Trials', num2str(i),'Fits.png'));
end
%end

%% Using only the data at 100 levels 
n_trial_this_cell_power = sum(powers_trials_this_cell==100);
responses_power = zeros(n_trial_this_cell_power, length(current_template));
stims_power = zeros(n_selected_cells,n_trial_this_cell_power, length(current_template));
loc_power=zeros(n_trial_this_cell_power,1);
power_power=zeros(n_trial_this_cell_power,1);
counter=0;

for i_trial = 1:n_trial_this_cell
    if powers_trials_this_cell(i_trial) == 100
        counter=counter+1;
    for i_cell = 1:n_selected_cells
        k=stim_this_cell(i_trial,i_cell);
        stim = current_template*k;
        stims_power(i_cell,counter,:) = stim;
        
    end
    loc_power(counter)=locations_trials_this_cell(i_trial,1);
    power_power(counter)=powers_trials_this_cell(i_trial);
    
    spikes_delay =round(mpp_this_cell(i_trial).times);
    responses_power(counter,spikes_delay)=1;
    end
end

%% Fit the model 
responses_reg=responses_power;
stims_reg=stims_power;
% Only assigh values in the first column 
gain_grid =[-0.01 0 0.05]+gain_crude;
initial_values = zeros(length(gain_grid),2*n_selected_cells);
initial_values(:,1)=gain_grid;
[stats_conv2,pred_prob2] = fit_lifglm_v6(responses_reg, stims_reg,in_params,...
    background_rate,v_reset_known,v_th_known,first_only,loss_type,...
    gain_lower_bound,gain_upper_bound,initial_values,delay_params,single_connection);
% Output:
% stats_conv2=stats_conv;
% pred_prob2=pred_prob;
%%
for i = 1:length(unique_locs)
    figure(i)
    mpp_this_loc = mpp_this_cell(locations_trials_this_cell(:,1)==unique_locs(i));
    powers_trials_this_loc = powers_trials_this_cell(locations_trials_this_cell(:,1)==unique_locs(i));
    plot(sort([mpp_this_loc(powers_trials_this_loc==100).times]),1:length([mpp_this_loc(powers_trials_this_loc==100).times]),...
        '.','col',colors(i),'MarkerSize',20)
    hold on;
    
    i_trial = min(find(loc_power(:,1) == unique_locs(i)) );
    prob_this_loc = pred_prob2(1,i_trial,:)+background_rate;
    prob_this_loc = reshape(prob_this_loc,[n_grid 1]);
    plot(1:n_grid,prob_this_loc*100,'col',colors(i))
    
    hold off;
    ylim([0,20]);
    xlim([0,300]);
    xlabel('Time (1/20 ms)');
    ylabel('Indices (irrelevant)');
    
    saveas(i,strcat('./Trials', num2str(i),'Fits100.png'));
end