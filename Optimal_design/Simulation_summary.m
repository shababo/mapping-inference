clear;
%%
N=200;
A_val = 200;
num_sim = 10;

dt_optimal_sum = zeros(200,num_sim);
NRE_conn_optimal_sum = zeros(N/10,num_sim);
NRE_amp_optimal_sum = zeros(N/10,num_sim);
AUC_conn_optimal_sum = zeros(N/10,num_sim);

dt_random_sum = zeros(200,num_sim);
NRE_conn_random_sum = zeros(N/10,num_sim);
NRE_amp_random_sum = zeros(N/10,num_sim);
AUC_conn_random_sum = zeros(N/10,num_sim);

dt_ortho_sum = zeros(200,num_sim);
NRE_conn_ortho_sum = zeros(N/10,num_sim);
NRE_amp_ortho_sum = zeros(N/10,num_sim);
AUC_conn_ortho_sum = zeros(N/10,num_sim);

outflnm = strcat('../Figures/Simulation-optimal/A', num2str(A_val)); 

for randomseed = 1:10
    
flnm=strcat('../Data/sim-results/A', num2str(A_val),'S', num2str(randomseed),'.mat'); 
    
load(flnm)

% Average computing time:


    dt_optimal_sum(:,randomseed) = time_record_optimal(2*(1:N)+1 ) - time_record_optimal(2*(1:N));
    dt_random_sum(:,randomseed) = time_record_random(2*(1:N)+1) - time_record_random(2*(1:N));
dt_ortho_sum(:,randomseed) = time_record_ortho(2*(1:N)+1) - time_record_ortho(2*(1:N));

% Average performance 

NRE_conn_optimal_sum(:,randomseed) = NRE_conn_optimal;
NRE_amp_optimal_sum(:,randomseed) = NRE_amp_optimal;
AUC_conn_optimal_sum(:,randomseed) = AUC_conn_optimal;



NRE_conn_random_sum(:,randomseed) = NRE_conn_random;
NRE_amp_random_sum(:,randomseed) = NRE_amp_random;
AUC_conn_random_sum(:,randomseed) = AUC_conn_random;


NRE_conn_ortho_sum(:,randomseed) = NRE_conn_ortho;
NRE_amp_ortho_sum(:,randomseed) = NRE_amp_ortho;
AUC_conn_ortho_sum(:,randomseed) = AUC_conn_ortho;
end 


%% Plotting 
figure(1)
plot(1:N, mean(dt_random_sum,2) ,'col',[1,0,0,1],'Linewidth',4);
hold on;
plot(1:N, mean(dt_optimal_sum,2) ,'col',[0,0,1,1],'Linewidth',4);
ylim([0,12]);
xlim([0,N]);
for i = 1:num_sim
plot(1:N, dt_random_sum(:,i),'col',[1,0,0,0.1],'Linewidth',1);
plot(1:N, dt_optimal_sum(:,i),'col',[0,0,1,0.1],'Linewidth',1);
end    

xlabel('Number of batches');
ylabel('Computing time per batch (seconds)');
hold off;

% NRE of connectivity reconstruction 

figure(2)
plot(1:(N/10), mean(NRE_conn_random_sum,2) ,'col',[1,0,0,1],'Linewidth',4);
hold on;
plot(1:(N/10), mean(NRE_conn_optimal_sum,2) ,'col',[0,0,1,1],'Linewidth',4);
xlim([0,N/10+1]);
ylim([0.3,2]);

for i = 1:num_sim
    plot(1:(N/10), NRE_conn_random_sum(:,i),'col',[1,0,0,0.1],'Linewidth',1);
    plot(1:(N/10), NRE_conn_optimal_sum(:,i),'col',[0,0,1,0.1],'Linewidth',1);
end    
hold off;

xlabel('Number of batches');
ylabel('NRE of connectivity');

xticks([0 5 10 15 20])
xticklabels({'0', '50', '100', '150', '200'})
hold off;



figure(3)
plot(1:(N/10), mean(NRE_amp_random_sum,2) ,'col',[1,0,0,1],'Linewidth',4);
hold on;
plot(1:(N/10), mean(NRE_amp_optimal_sum,2) ,'col',[0,0,1,1],'Linewidth',4);
xlim([0,N/10+1]);
ylim([0.3,2]);

for i = 1:num_sim
    plot(1:(N/10), NRE_amp_random_sum(:,i),'col',[1,0,0,0.1],'Linewidth',1);
    plot(1:(N/10), NRE_amp_optimal_sum(:,i),'col',[0,0,1,0.1],'Linewidth',1);
end    
hold off;

xlabel('Number of batches');
ylabel('NRE of amplitudes');

xticks([0 5 10 15 20])
xticklabels({'0', '50', '100', '150', '200'})
hold off;



figure(4)
plot(1:(N/10), mean(AUC_conn_random_sum,2) ,'col',[1,0,0,1],'Linewidth',4);
hold on;
plot(1:(N/10), mean(AUC_conn_optimal_sum,2) ,'col',[0,0,1,1],'Linewidth',4);
xlim([0,N/10+1]);
ylim([0.5,1]);

for i = 1:num_sim
    plot(1:(N/10), AUC_conn_random_sum(:,i),'col',[1,0,0,0.1],'Linewidth',1);
    plot(1:(N/10), AUC_conn_optimal_sum(:,i),'col',[0,0,1,0.1],'Linewidth',1);
end    
hold off;

xlabel('Number of batches');
ylabel('AUC of connectivity');
xticks([0 5 10 15 20])
xticklabels({'0', '50', '100', '150', '200'})
hold off;

saveas(1,strcat(outflnm,'Time','.jpg'));

saveas(2,strcat(outflnm,'NRE_conn','.jpg'));

saveas(3,strcat(outflnm,'NRE_amp','.jpg'));

saveas(4,strcat(outflnm,'AUC_conn','.jpg'));



%% Draw the plot for the ''precise'' stimulation


%% Visualize the stimulation locations
%% Parameters for the data generating mechanism
rng(12242,'twister');
%rng(12212,'twister');

% Parameters on the grid and experimental design 
region_width = 500;
region_height = 500;

% load parameters for the model
run('Parameters_setup_ground_truth.m')

% Pre-processing for  experiment
num_samples=50; % Number of samples to estimate the expected entropies
num_trials_batch=20;

num_peaks = 20;
% We consider new locations on a grid (rather than the cell somas) in order
% to gain extra information (e.g., to be able to distinguish two
% neighbouring cells)

num_dense=40; % Grid density

N=200; % Number of batches
%N_initial = 20;
num_trials_first = 200;

%forget_scale = 0.5;
sqrt_transform = false;
%design = 'random'; 
num_threshold=10;

A = diag([200, 200, 750]);
    num_sources = 4;  % Number of locations to stimulate in each trial
    run('Parameters_setup_experiment.m')
    % Run analysis and design

    
figure(11)
for i = 1:num_layers
    %connected_neurons_ind = find(neuron_features(i).amplitude);
    connected_neurons_ind = 1:size(neuron_features(i).amplitude,1);
    
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        (neuron_features(i).amplitude(connected_neurons_ind)+1)*35);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,0.3);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

xlim([20,460]);
ylim([-900,-400]);

loc_trials = reshape([location_random( (num_trials_first+1) :end,:) ], [(size(location_random,1)-num_trials_first)*4,1]);
loc_trials_counts  = tabulate(loc_trials);
loc_trials_counts_nonzero = loc_trials_counts(loc_trials_counts(:,2)>0,1:2);


stimulate = scatter(Z_dense(loc_trials_counts_nonzero(:,1),1),...
    -Z_dense(loc_trials_counts_nonzero(:,1),2),...
    loc_trials_counts_nonzero(:,2)*3, 'filled','d');
set(stimulate,'MarkerFaceColor','r');
alpha(stimulate,0.4);
hold on;

figure(12)

for i = 1:num_layers
    %connected_neurons_ind = find(neuron_features(i).amplitude);
    connected_neurons_ind = 1:size(neuron_features(i).amplitude,1);
    
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        (neuron_features(i).amplitude(connected_neurons_ind)+1)*35);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,0.3);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

xlim([20,460]);
ylim([-900,-400]);

loc_trials = reshape([location_optimal( (num_trials_first+1) :end,:) ], [(size(location_optimal,1)-num_trials_first)*4,1]);

loc_trials_counts  = tabulate(loc_trials);
loc_trials_counts_nonzero = loc_trials_counts(loc_trials_counts(:,2)>0,1:2);


stimulate = scatter(Z_dense(loc_trials_counts_nonzero(:,1),1),...
    -Z_dense(loc_trials_counts_nonzero(:,1),2),...
    loc_trials_counts_nonzero(:,2)*3, 'filled','d');
set(stimulate,'MarkerFaceColor','b');
alpha(stimulate,0.4);
hold on;
saveas(11,strcat(outflnm,'RandomDesign','.jpg'));
saveas(12,strcat(outflnm,'OptimalDesign','.jpg'));


% figure(13)
% 
% for i = 1:num_layers
%     %connected_neurons_ind = find(neuron_features(i).amplitude);
%     connected_neurons_ind = 1:size(neuron_features(i).amplitude,1);
%     
%     temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
%         -neuron_locations{i}(connected_neurons_ind,2),...
%         (neuron_features(i).amplitude(connected_neurons_ind)+1)*35);
%     set(temp,'MarkerFaceColor','k');
%     alpha(temp,0.3);
%     hold on
% end
% set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
% 
% xlim([20,460]);
% ylim([-900,-400]);
% 
% loc_trials = reshape([location_ortho( (num_trials_first+1) :end,:) ], [(size(location_ortho,1)-num_trials_first)*4,1]);
% loc_trials_counts  = tabulate(loc_trials);
% loc_trials_counts_nonzero = loc_trials_counts(loc_trials_counts(:,2)>0,1:2);
% 
% 
% stimulate = scatter(Z_dense(loc_trials_counts_nonzero(:,1),1),...
%     -Z_dense(loc_trials_counts_nonzero(:,1),2),...
%     loc_trials_counts_nonzero(:,2)*3, 'filled','d');
% set(stimulate,'MarkerFaceColor','g');
% alpha(stimulate,0.4);
% hold on;




%% Visualize the total effects on neurons 
XtX_optimal = X_optimal'*X_optimal;
eigs_optimal =eig(XtX_optimal);

XtX_random = X_random'*X_random;
eigs_random =eig(XtX_random);

XtX_ortho = X_ortho'*X_ortho;
eigs_ortho =eig(XtX_ortho);

figure(13)
%histogram(eigs_random)
histogram(diagsq(X_random))

figure(14)
%histogram(eigs_optimal)
histogram(diagsq(X_optimal))


figure(15)
%histogram(eigs_optimal)
histogram(diagsq(X_ortho))







%% Draw plots that include the orthogonal design 
%% Plotting 
figure(1)
plot(1:N, mean(dt_random_sum,2) ,'col',[1,0,0,1],'Linewidth',4);
hold on;
plot(1:N, mean(dt_optimal_sum,2) ,'col',[0,0,1,1],'Linewidth',4);
ylim([0,12]);
xlim([0,N]);
for i = 1:num_sim
plot(1:N, dt_random_sum(:,i),'col',[1,0,0,0.1],'Linewidth',1);
plot(1:N, dt_optimal_sum(:,i),'col',[0,0,1,0.1],'Linewidth',1);
end    

xlabel('Number of batches');
ylabel('Computing time per batch (seconds)');
hold off;

% NRE of connectivity reconstruction 

figure(2)
plot(1:(N/10), mean(NRE_conn_random_sum,2) ,'col',[1,0,0,1],'Linewidth',4);
hold on;
plot(1:(N/10), mean(NRE_conn_optimal_sum,2) ,'col',[0,0,1,1],'Linewidth',4);
plot(1:(N/10), mean(NRE_conn_ortho_sum,2) ,'col',[0,1,0,1],'Linewidth',4);

xlim([0,N/10+1]);
ylim([0.3,2]);

for i = 1:num_sim
    plot(1:(N/10), NRE_conn_random_sum(:,i),'col',[1,0,0,0.1],'Linewidth',1);
    plot(1:(N/10), NRE_conn_optimal_sum(:,i),'col',[0,0,1,0.1],'Linewidth',1);
    plot(1:(N/10), NRE_conn_ortho_sum(:,i),'col',[0,1,0,0.1],'Linewidth',1);
end    
hold off;

xlabel('Number of batches');
ylabel('NRE of connectivity');

xticks([0 5 10 15 20])
xticklabels({'0', '50', '100', '150', '200'})
hold off;



figure(3)
plot(1:(N/10), mean(NRE_amp_random_sum,2) ,'col',[1,0,0,1],'Linewidth',4);
hold on;
plot(1:(N/10), mean(NRE_amp_optimal_sum,2) ,'col',[0,0,1,1],'Linewidth',4);
plot(1:(N/10), mean(NRE_amp_ortho_sum,2) ,'col',[0,1,0,1],'Linewidth',4);

xlim([0,N/10+1]);
ylim([0.3,2]);

for i = 1:num_sim
    plot(1:(N/10), NRE_amp_random_sum(:,i),'col',[1,0,0,0.1],'Linewidth',1);
    plot(1:(N/10), NRE_amp_optimal_sum(:,i),'col',[0,0,1,0.1],'Linewidth',1);
    plot(1:(N/10), NRE_amp_ortho_sum(:,i),'col',[0,1,0,0.1],'Linewidth',1);
end    
hold off;

xlabel('Number of batches');
ylabel('NRE of amplitudes');

xticks([0 5 10 15 20])
xticklabels({'0', '50', '100', '150', '200'})
hold off;



figure(4)
plot(1:(N/10), mean(AUC_conn_random_sum,2) ,'col',[1,0,0,1],'Linewidth',4);
hold on;
plot(1:(N/10), mean(AUC_conn_optimal_sum,2) ,'col',[0,0,1,1],'Linewidth',4);
plot(1:(N/10), mean(AUC_conn_ortho_sum,2) ,'col',[0,1,0,1],'Linewidth',4);
xlim([0,N/10+1]);
ylim([0.5,1]);

for i = 1:num_sim
    plot(1:(N/10), AUC_conn_random_sum(:,i),'col',[1,0,0,0.1],'Linewidth',1);
    plot(1:(N/10), AUC_conn_optimal_sum(:,i),'col',[0,0,1,0.1],'Linewidth',1);
    plot(1:(N/10), AUC_conn_ortho_sum(:,i),'col',[0,1,0,0.1],'Linewidth',1);
end    
hold off;

xlabel('Number of batches');
ylabel('AUC of connectivity');
xticks([0 5 10 15 20])
xticklabels({'0', '50', '100', '150', '200'})
hold off;

%%
%% Visualization:
figure(90)


output_eva = output_random; 
run('Experiment_evaluate.m')
weighted_inclusion_random = overall_connectivity_w;


for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*35);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,0.3);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

xlim([20,460]);
ylim([-900,-400]);

overall_amplitudes(overall_amplitudes<0)=0;

    potential_neuron_grid = scatter(Z(:,1), -Z(:,2),...
        overall_amplitudes*35+0.1, 'filled','o');
    set(potential_neuron_grid,'MarkerFaceColor','r');
    alpha(potential_neuron_grid,0.3);
    hold on
    

% Minimize margin and remove ticks
axis off
set(gcf, 'Units','normal')
set(gca, 'Position',[0 0 1 1])
set(gca, 'xtick',[ ])
set(gca, 'ytick',[ ])
view(2)

hold off

saveas(90,strcat(outflnm,'amplitudes_map_random','.jpg'));

%%
figure(100)
output_eva = output_optimal; 
run('Experiment_evaluate.m')
weighted_inclusion_optimal= overall_connectivity_w;


for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*35);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,0.3);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

xlim([20,460]);
ylim([-900,-400]);

overall_amplitudes(overall_amplitudes<0)=0;

    potential_neuron_grid = scatter(Z(:,1), -Z(:,2),...
        overall_amplitudes*35+0.1, 'filled','o');
    set(potential_neuron_grid,'MarkerFaceColor','b');
    alpha(potential_neuron_grid,0.3);
    hold on
    

% Minimize margin and remove ticks
axis off
set(gcf, 'Units','normal')
set(gca, 'Position',[0 0 1 1])
set(gca, 'xtick',[ ])
set(gca, 'ytick',[ ])
view(2)

hold off

saveas(100,strcat(outflnm,'amplitudes_map_optimal','.jpg'));

