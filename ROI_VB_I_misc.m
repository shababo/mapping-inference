%% Visualizations for the analysis in Stage I

%% A map of all neurons 
figure(1)

for i = 1:num_layers
    temp = scatter(neuron_locations{i}(:,1),...
        -neuron_locations{i}(:,2),...
        25,'filled','o');
    set(temp,'MarkerFaceColor','k');
   alpha(temp,1);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

    xlim([20,460]);
    ylim([-900,-400]);

    temp = scatter(postsyn_position(:,1),...
        -postsyn_position(:,2),...
        80,'filled','o');
set(temp,'MarkerFaceColor','g');
   alpha(temp,1);

    % Minimize margin and remove ticks
  axis off
   set(gcf, 'Units','normal')
   set(gca, 'Position',[0 0 1 1])
   set(gca, 'xtick',[ ])
      set(gca, 'ytick',[ ])
      
hold off
%view(2)
saveas(1,'../Figures/Neurons_all.jpg')
%% A map of connected neurons 
figure(2)

for i = 1:num_layers
     connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*25,...
        'filled','o');
   
    set(temp,'MarkerFaceColor','k');
   alpha(temp,1);
    hold on
end
%set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

    xlim([20,460]);
    ylim([-900,-400]);

    temp = scatter(postsyn_position(:,1),...
        -postsyn_position(:,2),...
        182,'filled','o');
set(temp,'MarkerFaceColor','g');
   alpha(temp,1);

   
   % Minimize margin and remove ticks
   axis off
   set(gcf, 'Units','normal')
   set(gca, 'Position',[0 0 1 1])
   set(gca, 'xtick',[ ])
      set(gca, 'ytick',[ ])
hold off
%view(2)
saveas(2,'../Figures/Neurons_connected.jpg')

%---------------------%
%% A map of all neurons (shaded) and four locations 
figure(3)


     temp = scatter(Z(trial_locations_on_grid(1,:),1),...
        -Z(trial_locations_on_grid(1,:),2),...
        100,'filled','o');
set(temp,'MarkerFaceColor','b');
   alpha(temp,1);
hold on
for i = 1:num_layers
    
    temp = scatter(neuron_locations{i}(:,1),...
        -neuron_locations{i}(:,2),...
        25,'filled','o');
   
    set(temp,'MarkerFaceColor','k');
   alpha(temp,0.3);
    hold on
end
%set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

    xlim([20,460]);
    ylim([-900,-400]);

  
   
   % Minimize margin and remove ticks
   axis off
   set(gcf, 'Units','normal')
   set(gca, 'Position',[0 0 1 1])
   set(gca, 'xtick',[ ])
      set(gca, 'ytick',[ ])
hold off
%view(2)
saveas(3,'../Figures/Stage_I_trial.jpg')

%---------------------%

%% A map of pre-synaptic neurons and Stage I fit

figure(4)

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

    xlim([20,460]);
    ylim([-900,-400]);

for j = 1:size(amp_related_count_trials,2)
    coef = output(j).alpha;
    coef_thres = 0.7;
    potential_neuron_grid = scatter(Z(coef>coef_thres,1), -Z(coef>coef_thres,2), amplitude_threshold(j+1)*25,'filled','o');
    set(potential_neuron_grid,'MarkerFaceColor','r');
    alpha(potential_neuron_grid,0.4);
    hold on
    
end
   % Minimize margin and remove ticks
   axis off
   set(gcf, 'Units','normal')
   set(gca, 'Position',[0 0 1 1])
   set(gca, 'xtick',[ ])
      set(gca, 'ytick',[ ])
hold off
%saveas(4,'../Figures/Stage_I_fit.jpg')

%---------------------------------------------%
%% 

% params.A = A
params = struct;
params.A=A;
params.coords=Z(:,1:3);
params.K = size(Z,1);
params.N=N;

data=struct;
data.stims = trial_locations_on_grid;

% Unknows: 
params.eta = zeros(params.K,1);
params.sigma_s = ones(params.K,1);
params.sigma_n = 1;

params.t = 1:1:data_params.T;
params.tau = 10;
params.g = 1;
alpha_sum = sum(alpha_synapse(params.t,0,params.tau,-params.g));


pi_kr = exp(-0.5*squareform(pdist(params.coords,'mahalanobis',params.A)).^2);

pi_nk = zeros(params.N,params.K);
for n = 1:params.N
    pi_nk(n,:) = min(1,sum(pi_kr(:,data.stims(n,:)),2)');
end

outputall=struct([]);
    Y_n = sum(amp_related_count_trials,2);
    
    %Y_n = sum(data.responses,2)/alpha_sum;
    
    
    hyperparam_sigma_n = sqrt(length(params.t))*params.sigma_n/abs(alpha_sum);
    
    hyperparam_p_connected = .1*ones(params.K,1);
    
    alphas = zeros(params.K,1); %ones(params.K, 1) * alpha_0;
    mu = zeros(params.K, 1);
    s_sq = zeros(params.K,1);
    n_varbvs_samples = 5;
    % run_varbvs(X, Y, sigma_n, sigma_s, alpha, options)
    % run_varbvs_general(X, Y, sigma_n, sigma_s, alpha, eta, options);
    for sample = 1:n_varbvs_samples
        %[alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs(pi_nk>rand(params.N,params.K), Y_n, hyperparam_sigma_n, params.sigma_s(1), hyperparam_p_connected(1));%, params.eta);
        [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs_general(pi_nk>rand(params.N,params.K), Y_n, hyperparam_sigma_n, params.sigma_s(1), hyperparam_p_connected(1), params.eta);
        alphas = alphas+alpha_tmp/n_varbvs_samples;
        mu = mu+mu_tmp/n_varbvs_samples;
        s_sq = s_sq+s_sq_tmp/n_varbvs_samples;
    end
    
    outputall(1).alpha = alphas;
    outputall(1).mu = mu;
    outputall(1).s_sq = mu;
    
    %output(j).pi_kr = pi_kr;
    %output(j).pi_nk = pi_nk;
    % output(j).Y_scalar = Y_n;
    
    outputall(1).w_estimate = alphas.*mu;
%%
figure(90)

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

    coef = outputall.alpha;
    coef_thres = quantile(coef,0.80);
    potential_neuron_grid = scatter(Z(coef>coef_thres,1), -Z(coef>coef_thres,2), 45,'filled','o');
    set(potential_neuron_grid,'MarkerFaceColor','r');
    alpha(potential_neuron_grid,0.8);
    hold on

hold off
view(2)


%% Visualize sites to stimulate 
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
    alpha(temp,0.5);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

%selected_pixels = zeros(num_dense,num_dense,size(amp_related_count_trials,2));
%cent = cell(size(amp_related_count_trials,2),1);

    xlim([20,460]);
    ylim([-900,-400]);
    
     potential_neuron_grid = scatter(Z(:,1),...
     -Z(:,2),20,colormap(2,:),...
    'filled','d');
    set(potential_neuron_grid,'MarkerFaceColor','g');
    alpha(potential_neuron_grid,0.5);

    %x_range = ceil(sqrt(R)/(x_dense(2)-x_dense(1)));
%y_range = ceil(sqrt(R)/(y_dense(2)-y_dense(1)));

% for j = 1:size(amp_related_count_trials,2)
%        
%      potential_neuron_grid = scatter(next_sites_grid{j}(:,1),...
%      -next_sites_grid{j}(:,2),20,colormap(2,:),...
%     'filled','d');
%     set(potential_neuron_grid,'MarkerFaceColor','g');
%     alpha(potential_neuron_grid,0.5);
%      
%      potential_neuron_regional = scatter(next_sites_regional{j}(:,1),...
%      -next_sites_regional{j}(:,2),20,colormap(2,:),...
%     'filled','o');
%     set(potential_neuron_regional,'MarkerFaceColor','r');
%     alpha(potential_neuron_regional,0.5);
%      
%      potential_neuron_local = scatter(next_sites_local{j}(:,1),...
%          -next_sites_local{j}(:,2),amplitude_threshold(j)*10,colormap(2,:),...
%          'filled','o');
%      set(potential_neuron_local,'MarkerFaceColor','b');
%      alpha(potential_neuron_local,0.5);
% end

hold off
%saveas(10,'../Data/Sites_to_stimulate.jpg')
view(2)


%% Visualize the estimated coefficients:
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


%% Draw the inferred neurons and their synaptic strengths
% Merge the neurons in different layers into one big matrix:

all_likelihood_sum =zeros(0,size(amp_related_count_trials,2));
all_likelihood_max =zeros(0,size(amp_related_count_trials,2));
all_locations =zeros(0,3);
for i = 1:num_layers
    all_likelihood_sum = [all_likelihood_sum; neuron_likelihood_sum{i}];
    all_likelihood_max = [all_likelihood_max; neuron_likelihood_max{i}];
    all_locations = [all_locations; neuron_locations{i}(:,1:3)];
end

all_dist = diag(all_locations(:,1:2)*all_locations(:,1:2)')*ones(1,size(all_locations,1))...
    +ones(size(all_locations,1),1)*diag(all_locations(:,1:2)*all_locations(:,1:2)')'...
    -2*all_locations(:,1:2)*all_locations(:,1:2)';

%% Finding local maximum 
% We let the local maximum be the neurons with largests likelihood among
% its neighbours, given an acceptable radius R
% The distance matrix is calculated in the previous block 
R=2000;
figure(12)

colormap = jet(2);
colormap(1,:)= [1 1 1];
colormap(2,:)= [1 0 0];
for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*25);
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
    likelihood_thresholds =quantile( all_likelihood_sum(:,j), [0.80  0.9 0.99]);
    for i = 1:size(all_dist,1)
        if all_likelihood_sum(i,j) == max(all_likelihood_sum(all_dist(i,:)<R,j))
            if all_likelihood_sum(i,j)>likelihood_thresholds(3)
            neuron_centers = [neuron_centers; all_locations(i,:)];
            end
        end
    end
     potential_neuron = scatter(neuron_centers(:,1),-neuron_centers(:,2),amplitude_threshold(j+1)*25,colormap(2,:),'filled','o');
     alpha(potential_neuron,0.5);
end
hold off
saveas(12,'../Data/Final_I.jpg')

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
figure(1)
for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*50);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,.5);
    hold on
end

for j = 1:size(amp_related_count_trials,2)
    pvalues_grid = lmCount_related_amp{j}.Coefficients.pValue(1:end);
    %pvalues_grid =  lmCount_related_amp{j}.Coefficients.Estimate;
    
    a1=0.01;
    % P-values are proportional to the magnitudes
    %a1=quantile(pvalues_grid,0.2);
    %a2=quantile(pvalues_grid,0.3);
    %a1=quantile(pvalues_grid,0.95)
    
    xlim([20,460]);
    ylim([-900,-400]);
    %fullvote = scatter(Z(pvalues_grid>a1,1), -Z(pvalues_grid>a1,2),20,'filled','d');
    fullvote = scatter(Z(pvalues_grid<a1,1), -Z(pvalues_grid<a1,2),20,'filled','d');
    
    set(fullvote,'MarkerFaceColor','r');
    alpha(fullvote,1);
    
    hold on
    
    %saveas(j,'../Data/grid21j1.jpg')
end
spots = scatter(Z(:,1), -Z(:,2),10,'filled');
set(spots,'MarkerFaceColor','b');
alpha(spots,.2);

hold off
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
view(2)

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

% Visualize the Bayes estimates:

% figure(85)
% 
% for i = 1:num_layers
%     connected_neurons_ind = find(neuron_features(i).amplitude);
%     temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
%         -neuron_locations{i}(connected_neurons_ind,2),...
%         neuron_features(i).amplitude(connected_neurons_ind)*25);
%     set(temp,'MarkerFaceColor','k');
%    alpha(temp,0.8);
%     hold on
% end
% set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
% 
% %selected_pixels = zeros(num_dense,num_dense,size(amp_related_count_trials,2));
% %cent = cell(size(amp_related_count_trials,2),1);
% 
%     xlim([20,460]);
%     ylim([-900,-400]);
% %     
% %      potential_neuron_grid = scatter(Z(:,1),...
% %      -Z(:,2),20,colormap(2,:),...
% %     'filled','d');
% %     set(potential_neuron_grid,'MarkerFaceColor','k');
% %     alpha(potential_neuron_grid,0.2);
% 
% for j = 1:size(amp_related_count_trials,2)
%     coef = output(j).alpha;
%     coef_thres = 0.5;
%     potential_neuron_grid = scatter(Z(coef>coef_thres,1), -Z(coef>coef_thres,2), amplitude_threshold(j+1)*25,'filled','o');
%     set(potential_neuron_grid,'MarkerFaceColor','r');
%     alpha(potential_neuron_grid,0.4);
%     hold on
%     
% end
% 
% hold off
% view(2)
% 
% %% 
% % Visualize the estiamted weights 
% 
% % merging the estimated weights
% w_estimate_merge = [];
% mu_merge = [];
% alpha_merge = [];
% for j = 1:size(amp_related_count_trials,2)
%     w_estimate_merge = [w_estimate_merge output(j).w_estimate];
%     mu_merge = [mu_merge output(j).mu];
%     alpha_merge = [alpha_merge output(j).alpha];
% end
% 
% %% Obtain true amplitudes of the related neurons
% Z_amplitudes = all_amplitudes(neuron_in_region==1);
% 
% 
% 
% %% 
%    figure(11)
%    % find some cells that have non-zero amplitudes, and some have zero
%    % amplitudes
%    one_seq = 1:100;
%    nonzeros_seq = one_seq(Z_amplitudes(one_seq)>0);
%   chosen_ones = [1:5 nonzeros_seq(1:5)];
%    scale_y=0.5;
% for j = 1:length(chosen_ones)
%        
%     if Z_amplitudes(chosen_ones(j)) > 0
%         temp = line( amplitude_threshold(2:20),...
%             j/scale_y+2*w_estimate_merge(chosen_ones(j),:)/sum(w_estimate_merge(chosen_ones(j),:)),...
%             'LineWidth',3,...
%             'Color','k');
%         set(temp,'MarkerFaceColor','k');
%         %alpha(temp,0.8);
%         
%         %set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
%         
%         %selected_pixels = zeros(num_dense,num_dense,size(amp_related_count_trials,2));
%         %cent = cell(size(amp_related_count_trials,2),1);
%         
%         x = [amplitude_threshold(1):.1:amplitude_threshold(20)];
%         norm = normpdf(x,Z_amplitudes(chosen_ones(j)),evoked_params.sigma_a);
%         
%         line(x,...
%             j/scale_y+norm,...
%             'LineWidth',4,...
%             'Color','r');
%         xlim([amplitude_threshold(1),amplitude_threshold(20)]);
%         ylim([-0.1,length(chosen_ones)/scale_y+1.1]);
%         %
%     else 
%         if sum(w_estimate_merge(chosen_ones(j),:))/2<0.1
%             scale = 1;
%         else 
%             scale = sum(w_estimate_merge(chosen_ones(j),:));
%         end
%         
%         temp = line( amplitude_threshold(2:20),...
%             j/scale_y+2*w_estimate_merge(chosen_ones(j),:)/scale,...
%             'LineWidth',3,...
%             'Color','b');
%         %set(temp,'MarkerFaceColor','b');
%         %alpha(temp, 0.3);
%     end
% end
% hold off

%% Waste code 


% %% Fit regression by amplitudes 
% lmCount_related_amp=cell(size(amp_related_count_trials,2),1);
% %lmCount_related_amp_Robust=cell(size(amp_related_count_trials,2),1);
% for j = 1:size(amp_related_count_trials,2)
%     mdl_j=fitlm(covariates,amp_related_count_trials(:,j),'Intercept',false);
%     %EstCov = hac(mdl_j);
%     lmCount_related_amp{j}=mdl_j;
%     %lmCount_related_amp_Robust{j}=EstCov;
% end
% %% Visualize the estimates 
% % for j = 1:size(amp_related_count_trials,2)
% %     figure(j*17)
% %     histogram(lmCount_related_amp{j}.Coefficients.Estimate,30)
% %     
% %     % Check the quantiles and p-values 
% %     quantile(lmCount_related_amp{j}.Coefficients.Estimate,0.9)
% %     
% %     min(lmCount_related_amp{j}.Coefficients.Estimate(lmCount_related_amp{j}.Coefficients.pValue <0.05))
% %     
% % end
% figure(10)
% R=2500;
% 
% colormap = jet(2);
% colormap(1,:)= [1 1 1];
% colormap(2,:)= [1 0 0];
% for i = 1:num_layers
%     connected_neurons_ind = find(neuron_features(i).amplitude);
%     temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
%         -neuron_locations{i}(connected_neurons_ind,2),...
%         neuron_features(i).amplitude(connected_neurons_ind)*35);
%     set(temp,'MarkerFaceColor','k');
%     alpha(temp,0.8);
%     hold on
% end
% set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
% 
% selected_pixels = zeros(num_dense,num_dense,size(amp_related_count_trials,2));
% cent = cell(size(amp_related_count_trials,2),1);
% 
%     xlim([20,460]);
%     ylim([-900,-400]);
%     
%      potential_neuron_grid = scatter(Z(:,1),...
%      -Z(:,2),20,colormap(2,:),...
%     'filled','d');
%     set(potential_neuron_grid,'MarkerFaceColor','k');
%     alpha(potential_neuron_grid,0.2);
% 
% for j = 1:size(amp_related_count_trials,2)
%     coef = lmCount_related_amp{j}.Coefficients.Estimate;
%     coef_thres = 0.8;
%     potential_neuron_grid = scatter(Z(coef>coef_thres,1), -Z(coef>coef_thres,2), amplitude_threshold(j+1)*35,'filled','o');
%     set(potential_neuron_grid,'MarkerFaceColor','r');
%     alpha(potential_neuron_grid,0.4);
% hold on
%    
% end
% 
% hold off
% %saveas(10,'../Data/Sites_to_stimulate.jpg')
% view(2)
