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

% With intercept, it is rank deficient. 
covariates_intercept = [ones(num_combinations,1) covariates];
rank(covariates)
size(covariates)

%% Dividing the events by their amplitudes
% Now divide the events by quantiles of the amplitudes 
% We use overlapped regions to avoid separation due to 
num_threshold=20;
amplitude_threshold = quantile([related_mpp.amplitudes], (1/num_threshold)*[0:num_threshold]);
amp_related_count_trials = ones(size(trial_locations_on_grid,1),num_threshold-1);
for j = 1:(num_threshold-1)
    for i = 1:size(amp_related_count_trials,1)
        amp_related_count_trials(i,j) = sum(related_mpp(i).amplitudes>amplitude_threshold(j) & related_mpp(i).amplitudes<(amplitude_threshold(j+2)+0.01));
    end
end
%% Fit regression by amplitudes 
lmCount_related_amp=cell(size(amp_related_count_trials,2),1);
%lmCount_related_amp_Robust=cell(size(amp_related_count_trials,2),1);
for j = 1:size(amp_related_count_trials,2)
    mdl_j=fitlm(covariates,amp_related_count_trials(:,j),'Intercept',false);
    %EstCov = hac(mdl_j);
    lmCount_related_amp{j}=mdl_j;
    %lmCount_related_amp_Robust{j}=EstCov;
end
%% Visualize the estimates 
for j = 1:size(amp_related_count_trials,2)
    % figure(j*17)
    %histogram(lmCount_related_amp{j}.Coefficients.Estimate,30)
    
    % Check the quantiles and p-values 
    quantile(lmCount_related_amp{j}.Coefficients.Estimate,0.9)
    
    min(lmCount_related_amp{j}.Coefficients.Estimate(lmCount_related_amp{j}.Coefficients.pValue <0.05))
    
end

% Observations: 
% Using p-values as the criteria might pick sites that have negative
% estimates -- which corresponds to the unimportant sites...




%% Making some guestimates of where the neurons can be!
% Create a dense grid for evaluating the "likelihood" 
num_dense = 2e2+1;

likelihood = zeros(num_dense,num_dense,size(amp_related_count_trials,2));
x_dense = zeros(num_dense,1);
y_dense = zeros(num_dense,1);
x_dense = (0:(num_dense-1))*(max(Z(:,1))-min(Z(:,1)))/(num_dense-1) + min(Z(:,1));
y_dense = (0:(num_dense-1))*(max(Z(:,2))-min(Z(:,2)))/(num_dense-1) + min(Z(:,2));
Z_dense = zeros(num_dense^2,2);
% Take j=5 for example
for j = 1:size(amp_related_count_trials,2)
    j
    mdl_j=lmCount_related_amp{j};
    estimates = lmCount_related_amp{j}.Coefficients.Estimate;
    threshold = quantile(lmCount_related_amp{j}.Coefficients.Estimate,0.9);
    
    % we would need a better standard deviation estimates
    probability = mdl_j.Coefficients.Estimate(estimates > threshold);
    prob_sd = mdl_j.Coefficients.SE(estimates > threshold);
    locations = Z(estimates > threshold,:);
%     locations = Z_selected(pvalues_grid<0.05,:);
%     probability = mdl_j.Coefficients.Estimate;
%     prob_sd = mdl_j.Coefficients.SE;
%     locations = Z;
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

%% Find the regional maximum among the dense locations 
% There is a costumized threshold to be chosen...
threshold_level = 3;
quantile_level=0.95;
R=900;
x_range = ceil(sqrt(R)/(x_dense(2)-x_dense(1)));
y_range = ceil(sqrt(R)/(y_dense(2)-y_dense(1)));

next_sites_regional = cell(size(amp_related_count_trials,2),1);
next_sites_local = cell(size(amp_related_count_trials,2),1);
next_sites_grid = cell(size(amp_related_count_trials,2),1);

for j = 1:size(amp_related_count_trials,2)
    
    estimates = lmCount_related_amp{j}.Coefficients.Estimate;
    threshold = quantile(lmCount_related_amp{j}.Coefficients.Estimate,quantile_level);
    
    next_sites_grid{j} = Z(estimates > threshold,:);
    
    likelihood_thresholds =quantile( reshape(likelihood(:,:,j),[num_dense^2 1]), [0.90  0.95 0.99]);
    
    [cent{j}, varargout(:,:,j)]=FastPeakFind(likelihood(:,:,j),likelihood_thresholds(threshold_level));
    
    neuron_centers = zeros(size(cent{j},1)/2, 2);
    for i = 1:(size(cent{j},1)/2)
        neuron_centers(i,:) = [x_dense(cent{j}( 2*(i-1)+2)) y_dense(cent{j}( 2*(i-1)+1))]; 
    end
    
    threshold_region_vec= zeros(num_dense^2,1);
    for i = 1:num_dense
        for l = 1:num_dense
            threshold_region_vec((i-1)*num_dense + l) =varargout(i,l,j)+1;
        end
    end
    next_sites_local{j} = neuron_centers;
    neuron_centers=zeros(0,2);
    for i = 1:num_dense
        for l = 1:num_dense
            this_likelihood = likelihood(i,l,j);
            for ip = max(1, i-x_range) : min(num_dense, i+x_range)
                for lp = max(1, l-y_range) : min(num_dense, l+y_range)
                    dist = (x_dense(i)-x_dense(ip) )^2+ (y_dense(l)-y_dense(lp) )^2;
                        if dist < R
                            if this_likelihood < likelihood(ip,lp,j)
                                this_likelihood = likelihood(ip,lp,j);
                            end
                        end
                end
            end
            if this_likelihood < likelihood_thresholds(threshold_level)
                this_likelihood=0;
            end
            
            if this_likelihood == likelihood(i,l,j)
                neuron_centers = [neuron_centers; [x_dense(i) y_dense(l)] ];
            end
            
        end
    end
    next_sites_regional{j} = neuron_centers;
end

%% Run more experiments on these neurons
% On each site of interests, run ten trials.

for j = 1:size(amp_related_count_trials,2)
    size(next_sites_grid{j});
     size(next_sites_local{j});
       size(next_sites_regional{j});
end
merged_regional = cat(1,next_sites_regional{:});
size(merged_regional)
merged_regional = unique(merged_regional,'rows');
size(merged_regional)


merged_local = cat(1,next_sites_local{:});
size(merged_local)
merged_local = unique(merged_local,'rows');
size(merged_local)

merged_grid = cat(1,next_sites_grid{:});
size(merged_grid)
merged_grid = unique(merged_grid,'rows');
size(merged_grid)

%------------------------End of first stage-------------------------------------%
