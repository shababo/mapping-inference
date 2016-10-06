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

%% Loading functions and Data generation
addpath(genpath('../psc-detection'),genpath('../mapping-inference'),genpath('../mapping-core'));
%%
run('gendata_fullmodel_multi.m')


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
    auc_trials(i) = mean(Y_smooth);
end


denoised_auc_trials = ones(size(trial_locations_on_grid,1),1);
for i = 1:size(trial_locations_on_grid,1)
    denoised_auc_trials(i) = mean(mpp(i).trace);
end

count_trials = ones(size(trial_locations_on_grid,1),1);
for i = 1:size(trial_locations_on_grid,1)
    count_trials(i) = size(mpp(i).event_times,2);
end
% Use only the events in the first 400 time grid after the onsite of stimlus 
% evoked_params.stim_start
related_count_trials = ones(size(trial_locations_on_grid,1),1);
for i = 1:size(trial_locations_on_grid,1)
    related_count_trials(i) = sum(mpp(i).event_times>evoked_params.stim_start  & mpp(i).event_times< (400+evoked_params.stim_start));
end
covariates = zeros(size(trial_locations_on_grid,1), num_grids^2);


for i = 1:num_combinations
	covariates(i, trial_locations_on_grid(1 + num_repeats*(i-1),:)) = 1;    
end

covariates_intercept = [covariates ones(size(covariates,1),1)];

lmAUC=fitlm(covariates,auc_trials,'Intercept',false);

lmAUC_denoised=fitlm(covariates,denoised_auc_trials,'Intercept',false);

lmCount=fitlm(covariates,count_trials,'Intercept',false);

lmCount_related=fitlm(covariates,related_count_trials,'Intercept',false);





%% Visualizing the pvalues 
% Take the logorithm of the pvalues
  pvalues_grid = lmAUC.Coefficients.pValue(1:end);

a1=quantile(pvalues_grid,0.1);
a2=quantile(pvalues_grid,0.2);
figure(01)
for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        20*neuron_features(i).amplitude(connected_neurons_ind)*5,'.');
    hold on
end

fullvote = scatter(Z(pvalues_grid<a1,1), -Z(pvalues_grid<a1,2),80,'filled','d');
set(fullvote,'MarkerFaceColor','k');
alpha(fullvote,.4);

halfvote = scatter(Z(pvalues_grid<a2,1), -Z(pvalues_grid<a2,2),80,'filled','d');
set(halfvote,'MarkerFaceColor','b');
alpha(halfvote,.2);

%spots = scatter(Z(:,1), -Z(:,2),20,'filled');
%set(spots,'MarkerFaceColor','k');
%alpha(spots,.4);

hold off
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
view(2)


%% Now try the lasso regressions:
% Note that the default of glmnet is alpha = 1, which gives lasso penalty 
% addpath(genpath('../glmnet_matlab'));
% The glmnet crashes... try some examples to see why this happens


%% use the built-in method

glmnetAUC=lasso(covariates,auc_trials);

glmnetAUC_denoised=lasso(covariates,denoised_auc_trials);

glmnetCount=lasso(covariates,count_trials);

glmnetCount_related=lasso(covariates,related_count_trials);



%% Visualizing the solution path 
% Take the logorithm of the pvalues
solution_path = glmnetCount;
nonzeros = sum(solution_path~=0);

a1=max((nonzeros> (size(solution_path,1)*0.1)  ).*(1:100))+1;
a2=max((nonzeros> (size(solution_path,1)*0.2)  ).*(1:100))+1;
figure(10550)
for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        20*neuron_features(i).amplitude(connected_neurons_ind)*5,'.');
    hold on
end

fullvote = scatter(Z(solution_path(:,a1)~=0,1), -Z(solution_path(:,a1)~=0,2),80,'filled','d');
set(fullvote,'MarkerFaceColor','k');
alpha(fullvote,.4);

halfvote = scatter(Z(solution_path(:,a2)~=0,1), -Z(solution_path(:,a2)~=0,2),80,'filled','d');
set(halfvote,'MarkerFaceColor','b');
alpha(halfvote,.2);

%spots = scatter(Z(:,1), -Z(:,2),20,'filled');
%set(spots,'MarkerFaceColor','k');
%alpha(spots,.4);

hold off
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
view(2)


%% New thoughts:
temp=Y_grid{grid_locations(32,1),grid_locations(32,2)};

figure(12322)
num_stim = size(temp,1);
for i = 1:num_stim
    plot(1:2000,temp(i,:)+i*(20));
    hold on
end
hold off

%% 


temp=Y_grid{grid_locations(80,1),grid_locations(80,2)};

figure(12322)
num_stim = size(temp,1);
for i = 1:num_stim
    plot(1:2000,temp(i,:)+i*(20));
    hold on
end
hold off

%%
  pvalues_grid = lm1.Coefficients.pValue(2:end);


figure(12345)
for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        20*neuron_features(i).amplitude(connected_neurons_ind)*5,'.');
    hold on
end

fullvote = scatter(Z(32,1), -Z(32,2),80,'filled','d');
set(fullvote,'MarkerFaceColor','k');
alpha(fullvote,1);

halfvote = scatter(Z(80,1), -Z(80,2),80,'filled','d');
set(halfvote,'MarkerFaceColor','b');
alpha(halfvote,1);

hold off
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
view(2)


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
