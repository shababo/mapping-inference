%% Summary plots with raw data:
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Load the data set 
load('./Environments/6_3_s2c3_mpp_and_stim_data.mat')
target_locations =target_locs;
single_trial_limit=max(find(isnan(target_inds(:,2))));
trials_locations= target_inds(1:single_trial_limit,:);
trials_powers = stim_pow(1:single_trial_limit);
load('./Environments/6_3_cell_locs_reduced_and_fixed.mat')
cell_locations = cell_locs;
n_trial = size(trials_locations,1);
n_cell=size(cell_locations,1);
%%
max_time=300;
min_time=0;
for i_trial = 1:n_trial
    mpp(i_trial).locations = trials_locations(i_trial,:);
    mpp(i_trial).power = trials_powers(i_trial,:);

    if mpp(i_trial).num_events >0
        range_idx = mpp(i_trial).times<max_time & mpp(i_trial).times>min_time ;
        mpp(i_trial).num_events = sum(range_idx);
        mpp(i_trial).times = mpp(i_trial).times(range_idx);
        mpp(i_trial).amp =mpp(i_trial).amp(range_idx);
    end
end
%% Gibbs results
% Delay: 
%delay_params.type=2; %1: normal; 2: gamma
%delay_params.mean=60; delay_params.std=15;
load('./Results/Gibbs_samples_6_3.mat')%'gamma_samples','gain_samples','gamma_initial','gain_initial'

% gain_initial_temp=gain_samples(:,1);
% gamma_initial_temp=gamma_samples(:,1);


% gain_samples(:,1)=gain_initial_temp;
% gamma_samples(:,1)=gamma_initial_temp;

% Delay: 
%delay_params.type=2; %1: normal; 2: gamma
%delay_params.mean=60; delay_params.std=30;
% load('./Results/Gibbs_samples_6_3_v2.mat')%'gamma_samples','gain_samples','gamma_initial','gain_initial'

% 3:
%delay_params.type=2; %1: normal; 2: gamma
%delay_params.mean=35; delay_params.std=30;
% load('./Results/Gibbs_samples_6_3_v3.mat','gamma_samples','gain_samples','gamma_initial','gain_initial');

%%
gamma_fits = mean(gamma_samples,2);
gain_fits=mean(gain_samples,2);

% Count the number of events in trials at each locations (so that we can
% select the trials 
unique_locs_all =unique(trials_locations(:,1));
event_prop = zeros(length(unique_locs_all),1);
for i=1:length(unique_locs_all)
   relevant_index = trials_locations(:,1)==unique_locs_all(i); 
   event_prop(i) = length([mpp(relevant_index).times])/sum(relevant_index);
end

colors= ['w' 'g' 'b' 'k'];
alphas= [0.01 0.1 0.4 0.6];

prop_thresholds = [0 0.3 0.5];
colors_ind = zeros(length(unique_locs_all),1);
for i=1:length(unique_locs_all)
   colors_ind(i) = sum(event_prop(i)>prop_thresholds)+1;
end



figure(3)


for i=1:length(unique_locs_all)
temp3 = scatter(target_locs(unique_locs_all(i),2),-target_locs(unique_locs_all(i),1),...
    'Marker','s','SizeData',40,...
    'MarkerFaceColor','k', 'MarkerEdgeColor','w', 'MarkerFaceAlpha',alphas(colors_ind(i)));
hold on;
end

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])
% 
for i = 1:n_cell
    if gamma_fits(i)>0.05
    temp2 = scatter(cell_locations(i,2),-cell_locations(i,1),...
    400*gamma_fits(i),'MarkerEdgeColor','r','MarkerFaceColor','r',...
    'MarkerFaceAlpha',0.8*(cell_locations(i,3)/max(cell_locations(:,3))));
    end
hold on;
end

xlabel('X (um)');
ylabel('Y (um)');

axis off;
% Add legends to the plot 
temp5 = scatter([-150 -90 -30],[-155 -155 -155],...
    400*[0.8 0.5 0.2]);
set(temp5,'MarkerEdgeColor','r','MarkerFaceColor','r');
alpha(temp5,0.3);
txt4 = '\gamma = 0.8';
txt5 = '\gamma = 0.5';
txt6 = '\gamma = 0.2';
text(-140,-155,txt4)
text(-80,-155,txt5)
text(-20,-155,txt6)

xcord=[35 80 130];
ycord=[-155 -155 -155];
for i=1:3
temp5 = scatter(xcord(i),ycord(i),...
    'Marker','s','SizeData',40,...
     'MarkerFaceColor','k', 'MarkerEdgeColor','w', 'MarkerFaceAlpha',alphas(i+1));
end 
txt1 = 'N > 0';
txt2 = 'N > 0.3';
txt3 = 'N > 0.5';
text(40,-155,txt1)
text(85,-155,txt2)
text(135,-155,txt3)
% saveas(3,strcat('./Figures/Gibbs/','Gamma_final_v3','.jpg'));

%% Initial values v.s. final fits 

figure(4)


for i=1:length(unique_locs_all)
temp3 = scatter(target_locs(unique_locs_all(i),2),-target_locs(unique_locs_all(i),1),...
    'Marker','s','SizeData',40,...
    'MarkerFaceColor','k', 'MarkerEdgeColor','w', 'MarkerFaceAlpha',alphas(colors_ind(i)));
hold on;
end
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])
% 

for i = 1:n_cell
    if gamma_initial_temp(i)>0.05
    temp1 = scatter(cell_locations(i,2),-cell_locations(i,1),...
    400*gamma_initial_temp(i),'MarkerEdgeColor','b','MarkerFaceColor','b',...
    'MarkerFaceAlpha',0.8*(cell_locations(i,3)/max(cell_locations(:,3))));
    end
hold on;
end

xlabel('X (um)');
ylabel('Y (um)');

axis off;
% Add legends to the plot 
temp5 = scatter([-150 -90 -30],[-155 -155 -155],...
    400*[0.8 0.5 0.2]);
set(temp5,'MarkerEdgeColor','b','MarkerFaceColor','b');
alpha(temp5,0.3);
txt4 = '\gamma = 0.8';
txt5 = '\gamma = 0.5';
txt6 = '\gamma = 0.2';
text(-140,-155,txt4)
text(-80,-155,txt5)
text(-20,-155,txt6)

xcord=[35 80 130];
ycord=[-155 -155 -155];
for i=1:3
temp5 = scatter(xcord(i),ycord(i),...
    'Marker','s','SizeData',40,...
     'MarkerFaceColor','k', 'MarkerEdgeColor','w', 'MarkerFaceAlpha',alphas(i+1));
end 
txt1 = 'N > 0';
txt2 = 'N > 0.3';
txt3 = 'N > 0.5';
text(40,-155,txt1)
text(85,-155,txt2)
text(135,-155,txt3)
% saveas(4,strcat('./Figures/Gibbs/','Gamma_Initial','.jpg'));


%% Draw initial values against posterior means 
figure(5)

temp1 = scatter(gamma_initial_temp,gamma_fits,...
    30,'MarkerEdgeColor','k','MarkerFaceColor','k');
xlabel('Initial gammas');
ylabel('Posterior means');

xlim([min([gamma_initial_temp;gamma_fits] ) max([gamma_initial_temp;gamma_fits] )]);
ylim([min([gamma_initial_temp;gamma_fits] ) max([gamma_initial_temp;gamma_fits] )]);
saveas(5,strcat('./Figures/Gibbs/','Initial_vs_Final_gamma','.jpg'));


figure(6)
temp1 = scatter(gain_initial_temp,gain_fits,...
     30,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
xlabel('Initial gains');
ylabel('Posterior means');
xlim([min([gain_initial_temp;gain_fits] ) max([gain_initial_temp;gain_fits] )]);
ylim([min([gain_initial_temp;gain_fits] ) max([gain_initial_temp;gain_fits] )]);
saveas(6,strcat('./Figures/Gibbs/','Initial_vs_Final_gain','.jpg'));



%% Case study:

case_index=find((target_locs(:,2)>100 & target_locs(:,2)<130) & (target_locs(:,1)>30 & target_locs(:,1)<50));
case_index=case_index(1:2);


figure(4)


for i = 1:n_cell
    if gamma_fits(i)>0.05
    temp2 = scatter(cell_locations(i,2),-cell_locations(i,1),...
    400*gamma_fits(i),'MarkerEdgeColor','r','MarkerFaceColor','r',...
    'MarkerFaceAlpha',0.8*(cell_locations(i,3)/max(cell_locations(:,3))));
    end
hold on;
end

temp4 = scatter(target_locs(case_index,2),-target_locs(case_index,1),...
    'Marker','o','SizeData',100,...
    'MarkerFaceColor','g', 'MarkerEdgeColor','g','MarkerFaceAlpha',0.2);
hold on;

for i=1:length(unique_locs_all)
temp3 = scatter(target_locs(unique_locs_all(i),2),-target_locs(unique_locs_all(i),1),...
    'Marker','s','SizeData',40,...
    'MarkerFaceColor','k', 'MarkerEdgeColor','w', 'MarkerFaceAlpha',alphas(colors_ind(i)));
hold on;
end
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])

xlabel('X (um)');
ylabel('Y (um)');

axis off;
%  saveas(4,strcat('./Figures/Gibbs/','One_case','.jpg'));

%% Draw the responses
stim_counts = zeros(length(case_index),1);
relevant_trials_index = [];
for i_trial = 1:n_trial
    if sum(case_index== trials_locations(i_trial,1) )>0
        relevant_trials_index = [relevant_trials_index i_trial];
        stim_counts(find(case_index == trials_locations(i_trial,1)))=...
            stim_counts(find(case_index == trials_locations(i_trial,1)))+1;
    end
end

mpp_this_cell=mpp(relevant_trials_index);
trials_locations_this_cell = trials_locations(relevant_trials_index,:);
powers_trials_this_cell = trials_powers(relevant_trials_index);
length(relevant_trials_index)
%%
unique_locs= unique(trials_locations_this_cell(:,1));

colors=['r' 'g' 'b' 'k' 'y'];
for i = 1:length(unique_locs)
    figure(i)
    mpp_this_loc = mpp_this_cell(trials_locations_this_cell(:,1)==unique_locs(i));
    powers_trials_this_loc = powers_trials_this_cell(trials_locations_this_cell(:,1)==unique_locs(i));
    
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
    
    saveas(i,strcat('./Figures/Gibbs/Trials', num2str(unique_locs(i)),'.png'));
end
%end
%% Find the cell associated with this location:
cell_here=[];
for i = 1:length(unique_locs)
    dist=cell_locations-target_locations(unique_locs(i),:);
    sqdist=sum(dist.^2,2);
    cell_here = [cell_here find(sqdist<400)];
end
%% Case study: selected cells 
case_index=find((cell_locations(:,2)>50 & cell_locations(:,2)<80) & ...
    (cell_locations(:,1)>-20 & cell_locations(:,1)<0));

figure(4)

for i=1:length(unique_locs_all)
temp3 = scatter(target_locs(unique_locs_all(i),2),-target_locs(unique_locs_all(i),1),...
    'Marker','s','SizeData',40,...
    'MarkerFaceColor','k', 'MarkerEdgeColor','w', 'MarkerFaceAlpha',alphas(colors_ind(i)));
hold on;
end
for i = 1:n_cell
    if gamma_fits(i)>0.05
    temp2 = scatter(cell_locations(i,2),-cell_locations(i,1),...
    400*gamma_fits(i),'MarkerEdgeColor','r','MarkerFaceColor','r',...
    'MarkerFaceAlpha',0.8*(cell_locations(i,3)/max(cell_locations(:,3))));
    end
hold on;
end
% 
temp4 = scatter(cell_locations(case_index,2),-cell_locations(case_index,1),...
    'Marker','o','SizeData',200,...
    'MarkerFaceColor','g', 'MarkerEdgeColor','g','MarkerFaceAlpha',0.5);
hold on;



set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])
xlabel('X (um)');
ylabel('Y (um)');
axis off;
% saveas(4,strcat('./Figures/Gibbs/','Selected_cells','.jpg'));
%% Find targets near these locations 
target_here=[];
unique_locs=case_index;
for i = 1:length(unique_locs)
    dist=target_locations-cell_locations(unique_locs(i),:);
    sqdist=sum(dist.^2,2);
    target_here = [target_here find(sqdist<400)];
end

stim_counts = zeros(length(case_index),1);
relevant_trials_index = [];
for i_trial = 1:n_trial
    if sum(target_here== trials_locations(i_trial,1) )>0
        relevant_trials_index = [relevant_trials_index i_trial];
        stim_counts(find(target_here== trials_locations(i_trial,1)))=...
            stim_counts(find(target_here== trials_locations(i_trial,1)))+1;
    end
end

mpp_this_cell=mpp(relevant_trials_index);
trials_locations_this_cell = trials_locations(relevant_trials_index,:);
powers_trials_this_cell = trials_powers(relevant_trials_index);
length(relevant_trials_index)
%% Estimate the probability
gamma_fits(unique_locs)
gain_fits(unique_locs)
load('./Environments/chrome-template-3ms.mat');
downsamp=1;max_time=300;
power_level = unique(trials_powers);
num_power_level=length(power_level);
current_template=template(1:downsamp:max_time);
t_vect= 1:1:max_time;
% Load the shape template 
load('./Environments/l23_template_cell.mat');
%l23_average_shape
temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;
shape_template=l23_average_shape;

%%
delay_params.type=2; %1: normal; 2: gamma
delay_params.mean=30; delay_params.std=15;
delay_params.delayed=true; delay_params.n_grid=200;

delay_prob = zeros(delay_params.n_grid+1,1);
if delay_params.type == 1 % normal
    delay_prob = normpdf((0:delay_params.n_grid),delay_params.mean,...
        delay_params.std);
elseif delay_params.type == 2 % gamma
    shape=(delay_params.mean^2)/(delay_params.std^2);
    %scale
    scale = delay_params.mean/shape;
    delay_prob = gampdf((0:delay_params.n_grid),shape,scale);
end
% we approximate the probability with densities
delay_prob = delay_prob/sum(delay_prob);
min_delay = 0;
max_delay = delay_params.n_grid;

n_cell=size(cell_locations,1);
n_trial = length(mpp);
n_grid=length(current_template);

%%
%------------------------------%
% Calculate the size of stimuli
cell_params.locations =  cell_locations;
cell_params.shape_gain = ones(n_cell,1);
cell_template = struct();
cell_template.shape= shape_template;
[pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
    cell_template,target_locations);
stimuli_size=zeros(n_trial,n_cell);
for l = 1:n_trial
    for m = 1:size(mpp(l).locations,2)
        if isnan(mpp(l).locations(m))
        else
            stimuli_size(l,:) = stimuli_size(l,:)+...
                (pi_target(:,mpp(l).locations(m)).*mpp(l).power)';
        end
    end
end
%%
%------------------------------%
% Store the relevant trials for each cell
stim_threshold = 10;
g=0.02;
relevant_trials_per_cell=cell([n_cell 1]);
temp =1:n_trial;
for i_cell = 1:n_cell
    relevant_indicator=stimuli_size(:,i_cell)>stim_threshold;
    relevant_trials_per_cell{i_cell}=temp(relevant_indicator);
end

%------------------------------%
%   for each trial, save the sum log of probability of no spiking
%                        the sum probability of firing at each event
v_trace=cell([n_trial length(unique_locs)]);
for i_cell = 1:length(unique_locs)
    relevant_trials =relevant_trials_per_cell{unique_locs(i_cell)};
    for i = 1:length(relevant_trials)
        i_trial = relevant_trials(i);
        temp_trace = zeros([n_grid 1]);
        stims_temp=current_template*stimuli_size(i_trial,unique_locs(i_cell));
        temp_trace(1)=0;
        for i_t = 2:n_grid
            temp1=reshape(stims_temp(1:(i_t-1)), [i_t-1,1]);
            temp2=reshape(exp( ((1:(i_t-1)) -i_t+1)*g), [1,i_t-1]);
            temp_trace(i_t) = temp2*temp1;
        end
        v_trace{i_trial,i_cell}=temp_trace;
    end
    fprintf('Cell %d voltage grid done;\n', i_cell);
end

%%
v_th_known=15*ones(n_cell,1);
linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};
prob_sequences = cell(length(unique_locs), 1);
for i_cell = 1:length(unique_locs)
    relevant_trials =relevant_trials_per_cell{unique_locs(i_cell)};
    prob_sequence{i_cell}=cell(length(relevant_trials_per_cell{unique_locs(i_cell)}),1);
    for i = 1:length(relevant_trials)
        i_trial = relevant_trials(i);
        v_trace_one=v_trace{i_trial,i_cell}; gain_one=gain_fits(unique_locs(i_cell));
        v_th_known_one=v_th_known(unique_locs(i_cell));
        [prob_first_spike_delayed] = voltage_to_prob(gain_one,  v_trace_one,...
            v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob);
       prob_sequence{i_cell}{i_trial}=prob_first_spike_delayed;
    end
    fprintf('%d prob grid\n', i_cell);
end

%%

colors=['r' 'g' 'b' 'k' 'y'];
% for i = 1:length(target_here)
i=1;

    figure(i)
    trials_idx = find(trials_locations(:,1)==target_here(i));
    mpp_this_loc = mpp_this_cell(trials_locations_this_cell(:,1)==target_here(i));
    powers_trials_this_loc = powers_trials_this_cell(trials_locations_this_cell(:,1)==target_here(i));
    
    plot(sort([mpp_this_loc(powers_trials_this_loc==100).times]),1:length([mpp_this_loc(powers_trials_this_loc==100).times]),...
        '.','col','k','MarkerSize',20)
    hold on;
    plot(sort([mpp_this_loc(powers_trials_this_loc==75).times]),1:length([mpp_this_loc(powers_trials_this_loc==75).times]),...
        '.','col','k','MarkerSize',15)
    hold on;
    plot(sort([mpp_this_loc(powers_trials_this_loc==50).times]),1:length([mpp_this_loc(powers_trials_this_loc==50).times]),...
        '.','col','k','MarkerSize',10)
    
    
    i_trial = min(find(powers_trials_this_cell==75 ) );
    prob_this_loc = prob_sequence{1}{trials_idx(i_trial)};
    plot(1:n_grid,prob_this_loc*200,'col',colors(1),'LineStyle','--')
    hold on;
    i_trial = min(find(powers_trials_this_cell==100 ) );
    prob_this_loc = prob_sequence{1}{trials_idx(i_trial)};
    plot(1:n_grid,prob_this_loc*200,'col',colors(1))
    hold on;
    
    i_trial = min(find(powers_trials_this_cell==100 ) );
     prob_this_loc = prob_sequence{2}{trials_idx(i_trial)};
    plot(1:n_grid,prob_this_loc*200,'col',colors(2))
    hold on;
    i_trial = min(find(powers_trials_this_cell==75 ) );
     prob_this_loc = prob_sequence{2}{trials_idx(i_trial)};
    plot(1:n_grid,prob_this_loc*200,'col',colors(2),'LineStyle','--')
    
    
    hold off;
    ylim([0,20]);
    xlim([0,300]);
    xlabel('Time (1/20 ms)');
    ylabel('Indices (irrelevant)');
    
saveas(i,strcat('./Figures/Gibbs/Cells', num2str(unique_locs(i)),'.png'));

%% Parametric bootstrap 
load('./Results/Simulated_Gibbs_samples_6_3.mat')%'gamma_samples','gain_samples','gamma_initial','gain_initial'

%%
gamma_fits = mean(gamma_samples,2);
gain_fits=mean(gain_samples,2);

% Count the number of events in trials at each locations (so that we can
% select the trials 
unique_locs_all =unique(trials_locations(:,1));
event_prop = zeros(length(unique_locs_all),1);
for i=1:length(unique_locs_all)
   relevant_index = trials_locations(:,1)==unique_locs_all(i); 
   event_prop(i) = length([mpp(relevant_index).times])/sum(relevant_index);
end

colors= ['w' 'g' 'b' 'k'];
alphas= [0.01 0.1 0.4 0.6];

prop_thresholds = [0 0.3 0.5];
colors_ind = zeros(length(unique_locs_all),1);
for i=1:length(unique_locs_all)
   colors_ind(i) = sum(event_prop(i)>prop_thresholds)+1;
end

%%
figure(3)


for i=1:length(unique_locs_all)
temp3 = scatter(target_locs(unique_locs_all(i),2),-target_locs(unique_locs_all(i),1),...
    'Marker','s','SizeData',40,...
    'MarkerFaceColor','k', 'MarkerEdgeColor','w', 'MarkerFaceAlpha',alphas(colors_ind(i)));
hold on;
end

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])
% 
for i = 1:n_cell
    if gamma_fits(i)>0.05
    temp2 = scatter(cell_locations(i,2),-cell_locations(i,1),...
    400*gamma_fits(i),'MarkerEdgeColor','r','MarkerFaceColor','r',...
    'MarkerFaceAlpha',0.8*(cell_locations(i,3)/max(cell_locations(:,3))));
    end
hold on;
end

xlabel('X (um)');
ylabel('Y (um)');

axis off;
% Add legends to the plot 
temp5 = scatter([-150 -90 -30],[-155 -155 -155],...
    400*[0.8 0.5 0.2]);
set(temp5,'MarkerEdgeColor','r','MarkerFaceColor','r');
alpha(temp5,0.3);
txt4 = '\gamma = 0.8';
txt5 = '\gamma = 0.5';
txt6 = '\gamma = 0.2';
text(-140,-155,txt4)
text(-80,-155,txt5)
text(-20,-155,txt6)

xcord=[35 80 130];
ycord=[-155 -155 -155];
for i=1:3
temp5 = scatter(xcord(i),ycord(i),...
    'Marker','s','SizeData',40,...
     'MarkerFaceColor','k', 'MarkerEdgeColor','w', 'MarkerFaceAlpha',alphas(i+1));
end 
txt1 = 'N > 0';
txt2 = 'N > 0.3';
txt3 = 'N > 0.5';
text(40,-155,txt1)
text(85,-155,txt2)
text(135,-155,txt3)
saveas(3,strcat('./Figures/Gibbs/','Simulated_Gamma_final_v1','.jpg'));
%%

figure(4)


for i=1:length(unique_locs_all)
temp3 = scatter(target_locs(unique_locs_all(i),2),-target_locs(unique_locs_all(i),1),...
    'Marker','s','SizeData',40,...
    'MarkerFaceColor','k', 'MarkerEdgeColor','w', 'MarkerFaceAlpha',alphas(colors_ind(i)));
hold on;
end

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])
% 
for i = 1:n_cell
    if gamma_fits(i)>0.05
    temp2 = scatter(cell_locations(i,2),-cell_locations(i,1),...
    400*gamma_truth(i),'MarkerEdgeColor','r','MarkerFaceColor','r',...
    'MarkerFaceAlpha',0.8*(cell_locations(i,3)/max(cell_locations(:,3))));
    end
hold on;
end

xlabel('X (um)');
ylabel('Y (um)');

axis off;
% Add legends to the plot 
temp5 = scatter([-150 -90 -30],[-155 -155 -155],...
    400*[0.8 0.5 0.2]);
set(temp5,'MarkerEdgeColor','r','MarkerFaceColor','r');
alpha(temp5,0.3);
txt4 = '\gamma = 0.8';
txt5 = '\gamma = 0.5';
txt6 = '\gamma = 0.2';
text(-140,-155,txt4)
text(-80,-155,txt5)
text(-20,-155,txt6)

xcord=[35 80 130];
ycord=[-155 -155 -155];
for i=1:3
temp5 = scatter(xcord(i),ycord(i),...
    'Marker','s','SizeData',40,...
     'MarkerFaceColor','k', 'MarkerEdgeColor','w', 'MarkerFaceAlpha',alphas(i+1));
end 
txt1 = 'N > 0';
txt2 = 'N > 0.3';
txt3 = 'N > 0.5';
text(40,-155,txt1)
text(85,-155,txt2)
text(135,-155,txt3)
saveas(4,strcat('./Figures/Gibbs/','Simulated_Gamma_truth_v1','.jpg'));
%% Draw initial values against posterior means 
figure(5)

temp1 = scatter(gamma_truth,gamma_fits,...
    30,'MarkerEdgeColor','k','MarkerFaceColor','k');
xlabel('Initial gammas');
ylabel('Posterior means');

xlim([min([gamma_truth;gamma_fits] ) max([gamma_truth;gamma_fits] )]);
ylim([min([gamma_truth;gamma_fits] ) max([gamma_truth;gamma_fits] )]);
saveas(5,strcat('./Figures/Gibbs/','Simulated_Initial_vs_Final_gamma','.jpg'));


figure(6)
temp1 = scatter(gain_truth,gain_fits,...
     30,'MarkerEdgeColor','k','MarkerFaceColor','k');
hold on;
xlabel('Initial gains');
ylabel('Posterior means');
xlim([min([gain_truth;gain_fits] ) max([gain_truth;gain_fits] )]);
ylim([min([gain_truth;gain_fits] ) max([gain_truth;gain_fits] )]);
saveas(6,strcat('./Figures/Gibbs/','Simulated_Initial_vs_Final_gain','.jpg'));
