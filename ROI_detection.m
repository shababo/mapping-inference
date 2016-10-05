%% ROI detection: Stage I
% The goal of the first stage is to identify regions of interests
% This can yield a heatmap of possible synapses 

% Our strategy is to conduct a test for the num_repeats replicates on the
% same location. 
% The test can be one of these 
%   1) KS test against the homogeneous poisson process
%   2) Testing the mean (of time gaps) against the background level (p-value
%   =1 if no events)
%   3) Testing the mean counts against the background level
%   4) Mean PSC against background values 

%% Visualize the ground truth on a 2-d plane

figure(12345)
for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*5,'.');
    hold on
end

spots = scatter(Z(:,1), -Z(:,2),20,'filled');
set(spots,'MarkerFaceColor','k');
alpha(spots,.4);
hold off
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
view(2)
%% Data generation
run('gendata_fullmodel_multi.m')

%% KS test
% Obtain the background rate
count_trials = ones(size(trial_locations_on_grid,1),1);
for i = 1:size(trial_locations_on_grid,1)
    count_trials(i) = size(mpp(i).event_times,2);
end
% Two options:
muhat = mean(count_trials);
%muhat = median(count_trials)/log(2);

nulldist = makedist('Exponential','mu',2000/muhat/num_repeats); 

pvalues_grid = ones(num_combinations, num_grids^2);
for i = 1:num_combinations
    this_trial = trial_locations_on_grid(1 + num_repeats*(i-1),:);
    events_combined = mpp(1 + num_repeats*(i-1)).event_times;
    for nr = 2:num_repeats
        events_combined = [events_combined mpp(nr + num_repeats*(i-1)).event_times];
    end 
    events_combined = sort(events_combined);
    events_diff = events_combined(1);
    events_diff = [events_diff events_combined(2:end)-events_combined(1: end-1)];
    [h p k c] = kstest(events_diff,'CDF',nulldist);
    pvalues_grid(i, this_trial)= p;
end 

%% Tests on mean time gaps
% Obtain the background rate
% This is not working very well...
count_trials = ones(size(trial_locations_on_grid,1),1);
for i = 1:size(trial_locations_on_grid,1)
    count_trials(i) = size(mpp(i).event_times,2);
end
% Two options:
muhat = mean(count_trials);
%muhat = median(count_trials)/log(2);

% The mean gap is defined as the inverse of the muhat devided by the number
% of repeats
mean_gap = 1/muhat/num_repeats;

pvalues_grid = ones(num_combinations, num_grids^2);
for i = 1:num_combinations
    this_trial = trial_locations_on_grid(1 + num_repeats*(i-1),:);
    events_combined = mpp(1 + num_repeats*(i-1)).event_times;
    for nr = 2:num_repeats
        events_combined = [events_combined mpp(nr + num_repeats*(i-1)).event_times];
    end 
    events_combined = sort(events_combined);
    events_diff = events_combined(1);
    events_diff = [events_diff events_combined(2:end)-events_combined(1: end-1)];
    
    [h p] = ttest(events_diff,mean_gap);
    pvalues_grid(i, this_trial)= p;
end 

%% Tests on counts
% 
count_trials = ones(size(trial_locations_on_grid,1),1);
for i = 1:size(trial_locations_on_grid,1)
    count_trials(i) = size(mpp(i).event_times,2);
end
% Two options:
muhat = mean(count_trials);
%muhat = median(count_trials)/log(2);

pvalues_grid = ones(num_combinations, num_grids^2);
for i = 1:num_combinations
    this_trial = trial_locations_on_grid(1 + num_repeats*(i-1),:);
    count_this_trial = zeros(num_repeats,1); 
    for nr = 1:num_repeats
        count_this_trial(nr) = size(mpp(nr + num_repeats*(i-1)).event_times,2);
    end
    
    [h p] = ttest(count_this_trial,muhat);
    pvalues_grid(i, this_trial)= p;
end 


%% Tests on voltages
% Obtain the background Area Under the curve
% Might need to smooth the curves

auc_trials = ones(size(trial_locations_on_grid,1),1);
for i = 1:size(trial_locations_on_grid,1)
    Y_smooth = smooth(Y(i,:), 5,'moving');
    auc_trials(i) = mean(Y_smooth.^2);
end
% Two options:
muhat = mean(auc_trials);
%muhat = median(count_trials)/log(2);

pvalues_grid = ones(num_combinations, num_grids^2);
for i = 1:num_combinations
    this_trial = trial_locations_on_grid(1 + num_repeats*(i-1),:);
    auc_this_trial = auc_trials( (1:num_repeats) + num_repeats*(i-1));
    [h p] = ttest(auc_this_trial,muhat);
    pvalues_grid(i, this_trial)= p;
end 

%% Counting the votes
% Take the logorithm of the pvalues
log_pvalues_grid = log(pvalues_grid);
votes_on_grid = sum(pvalues_grid<0.05,1);
% Define colors based on the votes
% 0 = white 'w' [1 1 1]
% 1 = blue 'b' [0 0 1]
% 2 = black 'k' [0 0 0]
% color_spectrum = [1 1 1;1 1 1; 0 0 1; 0 0 0];
% colors_grid = color_spectrum(votes_on_grid+1,:);

figure(12345)
for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        10*neuron_features(i).amplitude(connected_neurons_ind)*5,'.');
    hold on
end

fullvote = scatter(Z(votes_on_grid>2,1), -Z(votes_on_grid>2,2),500,'filled');
set(fullvote,'MarkerFaceColor','k');
alpha(fullvote,.4);

halfvote = scatter(Z(votes_on_grid==2,1), -Z(votes_on_grid==2,2),500,'filled');
set(halfvote,'MarkerFaceColor','b');
alpha(halfvote,.2);
spots = scatter(Z(:,1), -Z(:,2),20,'filled');
set(spots,'MarkerFaceColor','k');
alpha(spots,.4);


hold off
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
view(2)

%% Regression
% The regression version takes one scalar response
% The covariate is a vector indicating which one is simulated 
% 
% We consider two responses here 
% - Counts over trials 
% - AUC of voltage 

auc_trials = ones(size(trial_locations_on_grid,1),1);
for i = 1:size(trial_locations_on_grid,1)
    Y_smooth = smooth(Y(i,:), 5,'moving');
    auc_trials(i) = mean(Y_smooth.^2);
end
% Two options:
muhat = mean(auc_trials);


count_trials = ones(size(trial_locations_on_grid,1),1);
for i = 1:size(trial_locations_on_grid,1)
    count_trials(i) = size(mpp(i).event_times,2);
end
muhat = mean(count_trials);

% Check the time to first event 
first_trials = ones(size(trial_locations_on_grid,1),1);
for i = 1:size(trial_locations_on_grid,1)
    if size(mpp(i).event_times,2)==0
        first_trials(i) = 2000;
    else 
        first_trials(i) = mpp(i).event_times(1);
    end
end

covariates = zeros(size(auc_trials,1), num_grids^2);
for i = 1:num_combinations
	covariates(i, trial_locations_on_grid(1 + num_repeats*(i-1),:)) = 1;    
end

covariates_intercept = [covariates ones(size(covariates,1),1)];

lm1=fitlm(covariates,auc_trials,'Intercept',false);

lm2=fitlm(covariates,count_trials,'Intercept',false);

lm3=fitlm(covariates,first_trials,'Intercept',false);


lm1.Coefficients.pValue(2:end);


%% Try an alternative way of doing KS test

% Obtain the background rate
count_trials = ones(size(trial_locations_on_grid,1),1);
for i = 1:size(trial_locations_on_grid,1)
    count_trials(i) = size(mpp(i).event_times,2);
end
% Two options:
muhat = mean(count_trials);
%muhat = median(count_trials)/log(2);

Y_grid = unstack_traces_multi(Y,trial_locations_on_grid, grid_locations);


mpp_grid = unstack_struct_multi(mpp,trial_locations_on_grid, grid_locations);

pvalues_grid = ones(num_grids^2,1);
count=1;
for i = 1:size(mpp_grid,1)
    for j = 1:size(mpp_grid,2)

        events_combined = mpp_grid{i,j}(1).event_times;
        for nr = 2:size(mpp_grid{i,j},2)
            events_combined = [events_combined mpp_grid{i,j}(nr).event_times];
        end 
        
        events_combined = sort(events_combined);
        events_diff = events_combined(1);
        events_diff = [events_diff events_combined(2:end)-events_combined(1: end-1)];
        
        nulldist = makedist('Exponential','mu',2000/muhat/size(mpp_grid{i,j},2)); 

        [h p k c] = kstest(events_diff,'CDF',nulldist);
        count = count +1;
        pvalues_grid(count)= p;
        p;
    end 
end


%% Visualizing the pvalues 
% Take the logorithm of the pvalues
  pvalues_grid = lm1.Coefficients.pValue(2:end);


figure(12345)
for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        20*neuron_features(i).amplitude(connected_neurons_ind)*5,'.');
    hold on
end

fullvote = scatter(Z(pvalues_grid<0.001,1), -Z(pvalues_grid<0.001,2),80,'filled','d');
set(fullvote,'MarkerFaceColor','k');
alpha(fullvote,.4);

halfvote = scatter(Z(pvalues_grid<0.01,1), -Z(pvalues_grid<0.01,2),80,'filled','d');
set(halfvote,'MarkerFaceColor','b');
alpha(halfvote,.2);

%spots = scatter(Z(:,1), -Z(:,2),20,'filled');
%set(spots,'MarkerFaceColor','k');
%alpha(spots,.4);

hold off
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
view(2)




