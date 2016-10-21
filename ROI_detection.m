%% ROI detection: Stage I
% The goal of the first stage is to identify regions of interests
% This can yield a heatmap of possible synapses 
%% Loading functions and Data generation
clear;
addpath(genpath('../psc-detection'),genpath('../mapping-inference'),genpath('../mapping-core'));
addpath(genpath('../Data'));

%%
run('gendata_fullmodel_multi.m')

%% Summary statistics to calculate 
% The regression version takes one scalar response
% The covariate is a vector indicating which one is simulated 
% 
% We consider two responses here 
% - Counts over trials 
% - AUC of voltage 

% auc_trials = ones(size(trial_locations_on_grid,1),1);
% for i = 1:size(trial_locations_on_grid,1)
%     Y_smooth = smooth(Y(i,:), 5,'moving');
%     auc_trials(i) = mean(Y_smooth);
% end
% 
% denoised_auc_trials = ones(size(trial_locations_on_grid,1),1);
% for i = 1:size(trial_locations_on_grid,1)
%     denoised_auc_trials(i) = mean(mpp(i).trace);
% end
% count_trials = ones(size(trial_locations_on_grid,1),1);
% for i = 1:size(trial_locations_on_grid,1)
%     count_trials(i) = size(mpp(i).event_times,2);
% end
% 
% amp_trials = ones(size(trial_locations_on_grid,1),1);
% for i = 1:size(trial_locations_on_grid,1)
%     amp_trials(i) = sum(mpp(i).amplitudes);
% end

% Use only the amplitudes in the related regions 
related_mpp = mpp;
unrelated_mpp = mpp;

for i = 1:size(trial_locations_on_grid,1)
    if size(mpp(i).event_times,2) > 0
        indices = mpp(i).event_times>evoked_params.stim_start  & mpp(i).event_times< (400+evoked_params.stim_start);
        related_mpp(i).amplitudes = mpp(i).amplitudes(indices);
        related_mpp(i).event_times = mpp(i).event_times(indices);
        unrelated_mpp(i).amplitudes = mpp(i).amplitudes(~indices);
        unrelated_mpp(i).event_times = mpp(i).event_times(~indices);
    end 
end

related_amp_trials = ones(size(trial_locations_on_grid,1),1);
for i = 1:size(trial_locations_on_grid,1)
    related_amp_trials(i) = sum(related_mpp(i).amplitudes);
end
% Use only the events in the first 400 time grid after the onsite of stimlus
related_count_trials = ones(size(trial_locations_on_grid,1),1);
for i = 1:size(trial_locations_on_grid,1)
    related_count_trials(i) = size(related_mpp(i).event_times,2);
end

unrelated_count_trials = ones(size(trial_locations_on_grid,1),1);
for i = 1:size(trial_locations_on_grid,1)
    unrelated_count_trials(i) = size(unrelated_mpp(i).event_times,2);
end

unrelated_amp_trials = ones(size(trial_locations_on_grid,1),1);
for i = 1:size(trial_locations_on_grid,1)
    unrelated_amp_trials(i) = sum(unrelated_mpp(i).amplitudes);
end

covariates = zeros(size(trial_locations_on_grid,1), num_grids^2);
for i = 1:num_combinations
	covariates(i, trial_locations_on_grid(1 + num_repeats*(i-1),:)) = 1;    
end


%% Run some regression 
lmCount_related=fitlm(covariates,related_count_trials,'Intercept',false);
% lmCount_unrelated=fitlm(covariates,unrelated_count_trials,'Intercept',false);
% lmAmp=fitlm(covariates,amp_trials,'Intercept',false);
% lmAmp_related=fitlm(covariates,related_amp_trials,'Intercept',false);
% lmAmp_unrelated=fitlm(covariates,unrelated_amp_trials,'Intercept',false);

% covariates_intercept = [covariates ones(size(covariates,1),1)];
% 
% lmAUC=fitlm(covariates,auc_trials,'Intercept',false);
% lmAUC_denoised=fitlm(covariates,denoised_auc_trials,'Intercept',false);
% 
% lmCount=fitlm(covariates,count_trials,'Intercept',false);



%% Dividing the events by their amplitudes
% Now divide the events by quantiles of the amplitudes 
num_threshold=15;
amplitude_threshold = quantile([related_mpp.amplitudes], (1/num_threshold)*[0:num_threshold]);
amp_related_count_trials = ones(size(trial_locations_on_grid,1),size(amplitude_threshold,2)-1);
for j = 1:size(amp_related_count_trials,2)
    for i = 1:size(amp_related_count_trials,1)
        amp_related_count_trials(i,j) = sum(related_mpp(i).amplitudes>amplitude_threshold(j) & related_mpp(i).amplitudes<(amplitude_threshold(j+1)+0.01));
    end
end
%% Run regression (linear/Poisson) on the category data 
% 
% Prepare the selected covariates 
pvalues_grid = lmCount_related.Coefficients.pValue(1:end);
selected_grid = [1:size(pvalues_grid,1)];
selected_grid = selected_grid(pvalues_grid<0.05);
covariates_selected = covariates(:,pvalues_grid<0.05);
Z_selected = Z(pvalues_grid<0.05,:);

lmCount_related_amp=cell(size(amp_related_count_trials,2),1);
for j = 1:size(amp_related_count_trials,2)
    mdl_j=fitlm(covariates_selected,amp_related_count_trials(:,j),'Intercept',false);
    lmCount_related_amp{j}=mdl_j;
end

%% Skip the first step, dive directly to the problem 
lmCount_related_amp=cell(size(amp_related_count_trials,2),1);
for j = 1:size(amp_related_count_trials,2)
    mdl_j=fitlm(covariates,amp_related_count_trials(:,j),'Intercept',false);
    lmCount_related_amp{j}=mdl_j;
end


%% Making some guestimates of where the neurons can be!
% Create a dense grid for evaluating the "likelihood" 
num_dense = 1e2+1;
likelihood = zeros(num_dense,num_dense,size(amp_related_count_trials,2));
x_dense = zeros(num_dense,1);
y_dense = zeros(num_dense,1);
x_dense = (0:(num_dense-1))*(max(Z(:,1))-min(Z(:,1)))/(num_dense-1) + min(Z(:,1));
y_dense = (0:(num_dense-1))*(max(Z(:,2))-min(Z(:,2)))/(num_dense-1) + min(Z(:,2));
Z_dense = zeros(num_dense^2,2);
% Take j=5 for example
for j = 1:size(amp_related_count_trials,2)
    mdl_j=lmCount_related_amp{j};
    pvalues_grid = lmCount_related_amp{j}.Coefficients.pValue(1:end);
    probability = mdl_j.Coefficients.Estimate(pvalues_grid<0.05);
    prob_sd = mdl_j.Coefficients.SE(pvalues_grid<0.05);
    locations = Z(pvalues_grid<0.05,:);
    % Draw a "heat map" based on the probability
    % We assume that the probability follows an Gaussian decay
    % The reference point is set to be the maximum of the estimated coefficient
    %
    reference_prob = min(max(probability),1);
    
    % Note here A is the scaling matrix!
    for i = 1:num_dense
        for l = 1:num_dense
            Z_dense((i-1)*num_dense + l,:) = [x_dense(i) y_dense(l)];
            for k = 1:size(probability,1)
                dist_scaled = (locations(k,1)-x_dense(i))^2/A(1,1)+(locations(k,2)- y_dense(l))^2/A(2,2);
                p_ijk = reference_prob*exp(-0.5*dist_scaled);
                % then check the Gaussian density
                likelihood(i,l,j) = likelihood(i,l,j)+normpdf(p_ijk,probability(k),prob_sd(k));
            end
        end
    end
end

%% Next, draw all estimated neurons in one plot
figure(10)

colormap = jet(2);
colormap(1,:)= [1 1 1];
colormap(2,:)= [1 0 0];
for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*10);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,0.5);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
view(2)


selected_pixels = zeros(num_dense,num_dense,size(amp_related_count_trials,2));
cent = cell(size(amp_related_count_trials,2),1);

    xlim([20,460]);
    ylim([-900,-400]);
    
    
for j = 1:size(amp_related_count_trials,2)
    likelihood_thresholds =quantile( reshape(likelihood(:,:,j),[num_dense^2 1]), [0.90  0.95 0.99]);
    [cent{j}, varargout(:,:,j)]=FastPeakFind(likelihood(:,:,j),likelihood_thresholds(2));
    
    neuron_centers = zeros(size(cent{j},1)/2, 2);
    for i = 1:(size(cent{j},1)/2)
        neuron_centers(i,:) = [y_dense(cent{j}( 2*(i-1)+1)) x_dense(cent{j}( 2*(i-1)+2))]; 
    end
    
    threshold_region_vec= zeros(num_dense^2,1);
    for i = 1:num_dense
        for l = 1:num_dense
            threshold_region_vec((i-1)*num_dense + l) =varargout(i,l,j)+1;
        end
    end
     potential_neuron = scatter(Z_dense(threshold_region_vec==2,1),-Z_dense(threshold_region_vec==2,2),amplitude_threshold(j)*10,colormap(2,:),'filled','o');
     alpha(potential_neuron,0.5);
end

hold off
%saveas(10,'../Data/guess41q95.jpg')

%% When the locations of neurons are known
% 
neuron_likelihood_sum = cell(num_layers,1);
neuron_likelihood_local_min = cell(num_layers,1);
for i = 1:num_layers
    neuron_likelihood_sum{i}= zeros(size(neuron_locations{i},1),size(amp_related_count_trials,2));
    neuron_likelihood_max{i}=zeros(size(neuron_locations{i},1),size(amp_related_count_trials,2));
end

num_layers;
neuron_locations{i}(:,1:2);
neuron_features(i).amplitude;

for j = 1:size(amp_related_count_trials,2)
    mdl_j=lmCount_related_amp{j};
    pvalues_grid = lmCount_related_amp{j}.Coefficients.pValue(1:end);
    probability = mdl_j.Coefficients.Estimate(pvalues_grid<0.05);
    prob_sd = mdl_j.Coefficients.SE(pvalues_grid<0.05);
    locations = Z(pvalues_grid<0.05,:);
    % Draw a "heat map" based on the probability
    % We assume that the probability follows an Gaussian decay
    % The reference point is set to be the maximum of the estimated coefficient
    %
    reference_prob = min(max(probability),1);
    
    % Note here A is the scaling matrix!
    for i = 1:num_layers
        for l = 1:size(neuron_locations{i},1)
            this_one = neuron_locations{i}(l,1:2);
            for k = 1:size(probability,1)
                dist_scaled = (locations(k,1)-this_one(1))^2/A(1,1)+(locations(k,2)- this_one(2))^2/A(2,2);
                p_ijk = reference_prob*exp(-0.5*dist_scaled);
                % then check the Gaussian density
                neuron_likelihood_sum{i}(l,j) = neuron_likelihood_sum{i}(l,j)+normpdf(p_ijk,probability(k),prob_sd(k));
                neuron_likelihood_max{i}(l,j) = max(neuron_likelihood_max{i}(l,j),normpdf(p_ijk,probability(k),prob_sd(k)));
                
            end
        end
    end
end

%% Draw the inferred neurons and their synaptic strengths
% Merge the neurons in different layers into one big matrix:
all_likelihood_sum =zeros(0,size(amp_related_count_trials,2));
all_likelihood_max =zeros(0,size(amp_related_count_trials,2));
all_locations =zeros(0,2);
for i = 1:num_layers
    all_likelihood_sum = [all_likelihood_sum; neuron_likelihood_sum{i}];
    all_likelihood_max = [all_likelihood_max; neuron_likelihood_max{i}];
    all_locations = [all_locations; neuron_locations{i}(:,1:2)];
end

all_dist = diag(all_locations*all_locations')*ones(1,size(all_locations,1))+ones(size(all_locations,1),1)*diag(all_locations*all_locations')'-2*all_locations*all_locations';

%% Finding local maximum 
% We let the local maximum be the neurons with largests likelihood among
% its neighbours, given an acceptable radius R
% The distance matrix is calculated in the previous block 
R=4900;
figure(11)

colormap = jet(2);
colormap(1,:)= [1 1 1];
colormap(2,:)= [1 0 0];
for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*10);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,0.6);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
view(2)

    xlim([20,460]);
    ylim([-900,-400]);
    
for j = 1:size(amp_related_count_trials,2)
    neuron_centers=zeros(0,2);
    likelihood_thresholds =quantile( all_likelihood_sum(:,j), [0.90  0.95 0.99]);
    for i = 1:size(all_dist,1)
        if all_likelihood_sum(i,j) == max(all_likelihood_sum(all_dist(i,:)<R,j))
            if all_likelihood_sum(i,j)>likelihood_thresholds(3)
            neuron_centers = [neuron_centers; all_locations(i,:)];
            end
        end
    end
     potential_neuron = scatter(neuron_centers(:,1),-neuron_centers(:,2),amplitude_threshold(j)*15,colormap(2,:),'filled','o');
     alpha(potential_neuron,0.5);
end
hold off
%saveas(10,'../Data/guess41q95.jpg')

%% Regression on distance matrix 
% Calculate the distance between the sites and the centers of neurons 
diff_mat=diag(all_locations*inv(A(1:2,1:2))*all_locations')*ones(1,size(Z,1))+(ones(size(all_locations,1),1)*diag(Z(:,1:2)*inv(A(1:2,1:2))*Z(:,1:2)')')-2*all_locations*inv(A(1:2,1:2))*Z(:,1:2)';

% Now calculate the design matrix:

covariates_neurons = zeros(size(trial_locations_on_grid,1), size(diff_mat,1));
for i = 1:size(trial_locations_on_grid,1)
	covariates_neurons(i,:) = sum(exp(-diff_mat(:,trial_locations_on_grid(1 + num_repeats*(i-1),:))/2),2);    
end
 
active_neurons = sum(covariates_neurons)>0.01;
active_locations = all_locations(active_neurons,:);
covariates_active_neurons = covariates_neurons(:,active_neurons);

lmCount_active_neurons=fitlm(covariates_active_neurons,related_count_trials,'Intercept',false);

%% Run subgroup regression
lmCount_active_amp=cell(size(amp_related_count_trials,2),1);
for j = 1:size(amp_related_count_trials,2)
    mdl_j=fitlm(covariates_active_neurons,amp_related_count_trials(:,j),'Intercept',false);
    lmCount_active_amp{j}=mdl_j;
end



%% Draw the crude p-values



pvalues_grid = lmCount_active_neurons.Coefficients.pValue(1:end);
%pvalues_grid =  -abs(lmCount_related.Coefficients.Estimate); 

% P-values are proportional to the magnitudes
%a1=quantile(pvalues_grid,0.2);
%a2=quantile(pvalues_grid,0.3);
a1=0.05;
a2=0.05;
figure(100)
for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*50);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,.5);
    hold on
end
xlim([20,460]);
ylim([-900,-400]);
fullvote = scatter(active_locations(pvalues_grid<a1,1), -active_locations(pvalues_grid<a1,2),20,'filled','d');
set(fullvote,'MarkerFaceColor','r');
alpha(fullvote,1);

halfvote = scatter(active_locations(pvalues_grid<a2,1), -active_locations(pvalues_grid<a2,2),20,'filled','d');
set(halfvote,'MarkerFaceColor','r');
alpha(halfvote,.5);


hold off
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
view(2)


%% Plotting subgroups

active_dist = diag(active_locations*active_locations')*ones(1,size(active_locations,1))+ones(size(active_locations,1),1)*diag(active_locations*active_locations')'-2*active_locations*active_locations';
R=2500;
    
    % P-values are proportional to the magnitudes
    %a1=quantile(pvalues_grid,0.2);
    %a2=quantile(pvalues_grid,0.3);
    a1=0.01;
   
    figure(101)
    for i = 1:num_layers
        connected_neurons_ind = find(neuron_features(i).amplitude);
        temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
            -neuron_locations{i}(connected_neurons_ind,2),...
            neuron_features(i).amplitude(connected_neurons_ind)*15);
        set(temp,'MarkerFaceColor','k');
        alpha(temp,.5);
        hold on
    end
    xlim([20,460]);
    ylim([-900,-400]);
for j = 1:size(amp_related_count_trials,2)
    neuron_centers=zeros(0,2);
    
    criteria=lmCount_active_amp{j}.Coefficients.pValue(1:end);
    for i = 1:size(active_dist,1)
        if criteria(i) == min(criteria(active_dist(i,:)<R))
            if criteria(i)<a1
            neuron_centers = [neuron_centers; active_locations(i,:)];
            end
        end
    end
     potential_neuron = scatter(neuron_centers(:,1),-neuron_centers(:,2),amplitude_threshold(j)*15,colormap(2,:),'filled','o');
     alpha(potential_neuron,0.5);
end

hold off
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
    view(2)
  

%% Data visualization

%% plot the neurons - colored by layer

%figure(123412)
%for i = 1:num_layers
%    scatter3(neuron_locations{i}(:,1),-neuron_locations{i}(:,2),neuron_locations{i}(:,3),'.');
%    hold on
%end
%hold off
%set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
%view(2)

% Draw the responses that cover this location 

%Y_grid = unstack_traces_multi(Y,trial_locations_on_grid, grid_locations);
% mpp_grid = unstack_struct(mpp,trial_grid_locations);
%plot_trace_stack_grid(Y_grid,10,1,0);

%%

% figure(12345)
%for i = 1:num_layers
%    connected_neurons_ind = find(neuron_features(i).amplitude);
%    scatter3(neuron_locations{i}(connected_neurons_ind,1),...
%        -neuron_locations{i}(connected_neurons_ind,2),...
%        neuron_locations{i}(connected_neurons_ind,3),...
%        neuron_features(i).amplitude(connected_neurons_ind)*5,'.');
%    hold on
%end
%hold off
%set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
%view(2)%%

% figure(123456)
% subplot(131)
% histogram(all_amplitudes(all_amplitudes ~= 0),1000)
% subplot(132)
% histogram(all_tau_rise(all_tau_rise ~= 0),1000)
% subplot(133)
% histogram(all_tau_fall(all_tau_fall ~= 0),1000)


%% Visualize the stimuli
% figure(12345)
% for i = 1:num_layers
%    connected_neurons_ind = find(neuron_features(i).amplitude);
%    scatter(neuron_locations{i}(connected_neurons_ind,1),...
%        -neuron_locations{i}(connected_neurons_ind,2),...
%        neuron_features(i).amplitude(connected_neurons_ind)*50,'.');
%    hold on
% end
% 
% spots = scatter(Z(:,1), -Z(:,2),20,'filled');
% set(spots,'MarkerFaceColor','k');
% alpha(spots,.4);
% hold off
% set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
% view(2)

%% Visualize the amplitudes by their period

figure(100)
all_amplitudes= [unrelated_mpp.amplitudes];
histogram(all_amplitudes,30,'Normalization','probability');
axis([0,35,0,0.2]);
title('Amplitudes of other events')

figure(111)
related_amplitudes= [related_mpp.amplitudes];
histogram(related_amplitudes,30,'Normalization','probability');
axis([0,35,0,0.2]);
title('Amplitudes of events between stimulus onset and stimulus onset + 20 ms')


%% Visualizing the pvalues 
% Take the logorithm of the pvalues
% Observations:
%   - Running linear regression using the related counts performs best
%       Reason is that the time period info is crucial 
%   - Running amplitudes is not very good 
%   - AUC performs similar, regardless of denoising 
pvalues_grid = lmCount_related.Coefficients.pValue(1:end);
%pvalues_grid =  -abs(lmCount_related.Coefficients.Estimate); 

% P-values are proportional to the magnitudes
%a1=quantile(pvalues_grid,0.2);
%a2=quantile(pvalues_grid,0.3);
a1=0.05;
a2=0.05;
figure(100)
for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*50);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,.5);
    hold on
end
xlim([20,460]);
ylim([-900,-400]);
fullvote = scatter(Z(pvalues_grid<a1,1), -Z(pvalues_grid<a1,2),20,'filled','d');
set(fullvote,'MarkerFaceColor','r');
alpha(fullvote,1);

halfvote = scatter(Z(pvalues_grid<a2,1), -Z(pvalues_grid<a2,2),20,'filled','d');
set(halfvote,'MarkerFaceColor','r');
alpha(halfvote,.5);

spots = scatter(Z(:,1), -Z(:,2),10,'filled');
set(spots,'MarkerFaceColor','b');
alpha(spots,.2);

hold off
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
view(2)
%saveas(1,'../Data/grid31.jpg')

%% Plotting subgroups
for j = 1:1
    pvalues_grid = lmCount_related_amp{j}.Coefficients.pValue(1:end);
    %pvalues_grid =  -abs(lmCount_related.Coefficients.Estimate);
    
    % P-values are proportional to the magnitudes
    %a1=quantile(pvalues_grid,0.2);
    %a2=quantile(pvalues_grid,0.3);
    a1=0.05;
    a2=0.05;
    figure(j)
    for i = 1:num_layers
        connected_neurons_ind = find(neuron_features(i).amplitude);
        temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
            -neuron_locations{i}(connected_neurons_ind,2),...
            neuron_features(i).amplitude(connected_neurons_ind)*50);
        set(temp,'MarkerFaceColor','k');
        alpha(temp,.5);
        hold on
    end
    xlim([20,460]);
    ylim([-900,-400]);
    fullvote = scatter(Z_selected(pvalues_grid<a1,1), -Z_selected(pvalues_grid<a1,2),20,'filled','d');
    set(fullvote,'MarkerFaceColor','r');
    alpha(fullvote,1);
    
    halfvote = scatter(Z_selected(pvalues_grid<a2,1), -Z_selected(pvalues_grid<a2,2),20,'filled','d');
    set(halfvote,'MarkerFaceColor','r');
    alpha(halfvote,.5);
    
    spots = scatter(Z(:,1), -Z(:,2),10,'filled');
    set(spots,'MarkerFaceColor','b');
    alpha(spots,.2);
    
    hold off
    set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
    view(2)
    %saveas(j,'../Data/grid21j1.jpg')
end

%% Visualizing 
colormap = jet(4);
colormap(1,:)= [1 1 1];
colormap(2,:)= [1 0.4 0.6];
colormap(3,:)= [1 0.4 0.6];
colormap(4,:)= [1 0 0];

for j = 1:1
    likelihood_thresholds =quantile( reshape(likelihood(:,:,j),[num_dense^2 1]), [0.8  0.95 0.99]);
    threshold_region_vec= zeros(num_dense^2,1);
    for i = 1:num_dense
        for l = 1:num_dense
            threshold_region_vec((i-1)*num_dense + l) =sum(likelihood(i,l,j)>likelihood_thresholds)+1;
        end
    end
    
    figure(j+100)
    potential_neuron = scatter(Z_dense(:,1),-Z_dense(:,2),20,colormap(threshold_region_vec,:),'filled','o');
    alpha(potential_neuron,1);
    
    xlim([20,460]);
    ylim([-900,-400]);
    
    hold on
    for i = 1:num_layers
        connected_neurons_ind = find(neuron_features(i).amplitude);
        temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
            -neuron_locations{i}(connected_neurons_ind,2),...
            neuron_features(i).amplitude(connected_neurons_ind)*5);
        set(temp,'MarkerFaceColor','k');
        alpha(temp,1);
        hold on
    end
    hold off
    set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
    view(2)
    saveas(j+100,'../Data/heat21j1.jpg')
end
% jet(10)
colormap = jet(4);
colormap(1,:)= [1 1 1];
colormap(2,:)= [1 1 1];
colormap(3,:)= [1 0.4 0.6];
colormap(4,:)= [1 0 0];

for j = 1:1
    likelihood_thresholds =quantile( reshape(likelihood(:,:,j),[num_dense^2 1]), [0.9  0.95 0.99]);
    threshold_region_vec= zeros(num_dense^2,1);
    for i = 1:num_dense
        for l = 1:num_dense
            threshold_region_vec((i-1)*num_dense + l) =sum(likelihood(i,l,j)>likelihood_thresholds)+1;
        end
    end
    
    figure(j+100)
    potential_neuron = scatter(Z_dense(:,1),-Z_dense(:,2),20,colormap(threshold_region_vec,:),'filled','o');
    alpha(potential_neuron,1);
    
    xlim([20,460]);
    ylim([-900,-400]);
    
    hold on
    for i = 1:num_layers
        connected_neurons_ind = find(neuron_features(i).amplitude);
        temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
            -neuron_locations{i}(connected_neurons_ind,2),...
            neuron_features(i).amplitude(connected_neurons_ind)*5);
        set(temp,'MarkerFaceColor','k');
        alpha(temp,1);
        hold on
    end
    hold off
    set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
    view(2)
    saveas(j+100,'../Data/thres21j1.jpg')
end
%%
colormap = jet(2);
colormap(1,:)= [1 1 1];
colormap(2,:)= [1 0 0];
selected_pixels = zeros(num_dense,num_dense,size(amp_related_count_trials,2));
cent = cell(size(amp_related_count_trials,2),1);
for j = 1:1 %size(amp_related_count_trials,2)
    likelihood_thresholds =quantile( reshape(likelihood(:,:,j),[num_dense^2 1]), [0.9  0.95 0.99]);
    
    [cent{j}, varargout(:,:,j)]=FastPeakFind(likelihood(:,:,j),likelihood_thresholds(2));
    
    neuron_centers = zeros(size(cent,1)/2, 2);
    for i = 1:(size(cent{j},1)/2)
        neuron_centers(i,:) = [y_dense(cent{j}( 2*(i-1)+1)) x_dense(cent{j}( 2*(i-1)+2))]; 
    end
    
    threshold_region_vec= zeros(num_dense^2,1);
    for i = 1:num_dense
        for l = 1:num_dense
            threshold_region_vec((i-1)*num_dense + l) =varargout(i,l,j)+1;
        end
    end
    
    
    figure(j)
     potential_neuron = scatter(Z_dense(threshold_region_vec==2,1),-Z_dense(threshold_region_vec==2,2),30,colormap(2,:),'filled','o');
     alpha(potential_neuron,1);
    
    xlim([20,460]);
    ylim([-900,-400]);
    
    hold on
    for i = 1:num_layers
        connected_neurons_ind = find(neuron_features(i).amplitude);
        temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
            -neuron_locations{i}(connected_neurons_ind,2),...
            neuron_features(i).amplitude(connected_neurons_ind)*5);
        set(temp,'MarkerFaceColor','k');
        alpha(temp,1);
        hold on
    end
    hold off
    set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
    view(2)
    saveas(j,'../Data/guess21j1.jpg')
end

