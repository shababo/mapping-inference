%% Script for KS test on single-location stimuli 
%% Loading functions and Data generation
clear;
addpath(genpath('../psc-detection'),genpath('../mapping-inference'),genpath('../mapping-core'));

%%
run('gendata_fullmodel.m')

%% KS tests

% Obtain the background rate from the sponetaneous events
count_trials = ones(size(trial_grid_locations,1),1);
for i = 1:size(trial_grid_locations,1)
    count_trials(i) = size(mpp(i).event_times,2);
end
% Two options, mean or median:
muhat = mean(count_trials);
%muhat = median(count_trials)/log(2);

% The rate is then 
muhat = muhat/(2000);


% Obtain the null distribution of the time gaps among stacked events
nulldist = makedist('Exponential','mu',1/muhat/num_repeats); 

num_grids=21;
pvalues_grid = ones(num_grids^2,1);
for i = 1:(num_grids^2)
    this_trial = trial_grid_locations(1 + num_repeats*(i-1),:);
    events_combined = mpp(1 + num_repeats*(i-1)).event_times;
    for nr = 2:num_repeats
        events_combined = [events_combined mpp(nr + num_repeats*(i-1)).event_times];
    end
    if size(events_combined,2)>0
        events_combined = sort(events_combined);
        events_diff = events_combined(1);
        events_diff = [events_diff events_combined(2:end)-events_combined(1: end-1)];
        [h p k c] = kstest(events_diff,'CDF',nulldist);
        pvalues_grid(this_trial(1),this_trial(2))= p;
        %pvalues_grid(i)= p;
    end
end 

%% Visualizing results
sign_index = pvalues_grid<0.05;
figure(12345)
for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        10*neuron_features(i).amplitude(connected_neurons_ind)*5,'.');
    hold on
end

for i = 1:num_grids
    for j = 1:num_grids
        if sign_index(i,j) == 1
            scatter((i-1)*20 - 200 + postsyn_position(1), -((j-1)*20 - 200 + postsyn_position(2)),50,'filled','k');
        end
    end
end


hold off
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
view(2)

