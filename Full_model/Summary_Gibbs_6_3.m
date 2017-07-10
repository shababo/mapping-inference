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

gamma_fits = mean(gamma_samples,2);
gain_samples=mean(gain_samples,2);

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



figure(1)


for i=1:length(unique_locs_all)
temp3 = scatter(target_locs(unique_locs_all(i),2),-target_locs(unique_locs_all(i),1),...
    'Marker','s','SizeData',40,...
    'MarkerFaceColor','k', 'MarkerEdgeColor','w', 'MarkerFaceAlpha',alphas(colors_ind(i)));
hold on;
end

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])

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
