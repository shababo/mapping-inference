%% ROI detection: Stage I
% If we can stimulate directly on the neurons
% 
%% Loading functions and Data generation
clear;
addpath(genpath('../psc-detection'),genpath('../mapping-inference'),genpath('../mapping-core'));
addpath(genpath('../Data'));
%%
run('gendata_fullmodel_multicells.m')
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

covariates = zeros(size(trial_locations_on_grid,1), size(Z,1));
for i = 1:N
	covariates(i, trial_locations_on_grid(i,:)) = 1;    
end

% With intercept, it is rank deficient. 
covariates_intercept = [ones(N,1) covariates];
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
% for j = 1:size(amp_related_count_trials,2)
%     figure(j*17)
%     histogram(lmCount_related_amp{j}.Coefficients.Estimate,30)
%     
%     % Check the quantiles and p-values 
%     quantile(lmCount_related_amp{j}.Coefficients.Estimate,0.9)
%     
%     min(lmCount_related_amp{j}.Coefficients.Estimate(lmCount_related_amp{j}.Coefficients.pValue <0.05))
%     
% end
figure(10)
R=2500;

colormap = jet(2);
colormap(1,:)= [1 1 1];
colormap(2,:)= [1 0 0];
for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*25);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,0.8);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

%selected_pixels = zeros(num_dense,num_dense,size(amp_related_count_trials,2));
%cent = cell(size(amp_related_count_trials,2),1);

    xlim([20,460]);
    ylim([-900,-400]);
%     
%      potential_neuron_grid = scatter(Z(:,1),...
%      -Z(:,2),20,colormap(2,:),...
%     'filled','d');
%     set(potential_neuron_grid,'MarkerFaceColor','k');
%     alpha(potential_neuron_grid,0.2);

for j = 1:size(amp_related_count_trials,2)
    coef = lmCount_related_amp{j}.Coefficients.Estimate;
    coef_thres = quantile(coef,0.98);
    potential_neuron_grid = scatter(Z(coef>coef_thres,1), -Z(coef>coef_thres,2), amplitude_threshold(j+1)*25,'filled','o');
    set(potential_neuron_grid,'MarkerFaceColor','r');
    alpha(potential_neuron_grid,0.4);
hold on
   
end

hold off
%saveas(10,'../Data/Sites_to_stimulate.jpg')
view(2)

%% The VB inference

% params.A = A
params = struct;
params.A=A;
params.coords=all_locations(:,1:3);
params.K = size(all_locations,1);
params.N=N;

data=struct;
data.stims = trial_locations_on_grid;

% Unknows: 
params.eta
params.sigma_s
params.sigma_n

params.t
params.tau
params.g

alpha_sum = sum(alpha_synapse(params.t,0,params.tau,-params.g));
hyperparam_sigma_n = sqrt(length(params.t))*params.sigma_n/abs(alpha_sum);


pi_kr = exp(-0.5*squareform(pdist(params.coords,'mahalanobis',params.A)).^2);

pi_nk = zeros(params.N,params.K);
for n = 1:params.N
    pi_nk(n,:) = min(1,sum(pi_kr(:,data.stims(n,:)),2)');
end
hyperparam_p_connected = .1*ones(params.K,1);




Y_n = sum(data.responses,2)/alpha_sum;

alpha = zeros(params.K,1); %ones(params.K, 1) * alpha_0;
mu = zeros(params.K, 1);
s_sq = zeros(params.K,1);
n_varbvs_samples = 20;
% run_varbvs(X, Y, sigma_n, sigma_s, alpha, options)
% run_varbvs_general(X, Y, sigma_n, sigma_s, alpha, eta, options);
for sample = 1:n_varbvs_samples
    %[alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs(pi_nk>rand(params.N,params.K), Y_n, hyperparam_sigma_n, params.sigma_s(1), hyperparam_p_connected(1));%, params.eta);
    [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs_general(pi_nk>rand(params.N,params.K), Y_n, hyperparam_sigma_n, params.sigma_s(1), hyperparam_p_connected(1), params.eta);
    alpha = alpha+alpha_tmp/n_varbvs_samples;
    mu = mu+mu_tmp/n_varbvs_samples;
    s_sq = s_sq+s_sq_tmp/n_varbvs_samples;
end

output.alpha = alpha;
output.mu = mu;
output.s_sq = mu;

output.pi_kr = pi_kr;
output.pi_nk = pi_nk;
output.Y_scalar = Y_n;


output.w_estimate = alpha.*mu;




%------------------------End of first stage-------------------------------------%
