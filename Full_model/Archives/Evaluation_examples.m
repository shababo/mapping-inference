%% Evaluation via simulations
%% Loading functions and Data generation
clear;
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Parameters for the data generating mechanism
rng(12242,'twister');

n_cell = 2;
n_trial = 1;

% load parameters for the model
% run('./Parameters_setup_ground_truth.m')
run('./Parameters_setup_single_LIF.m')
all_V_th = [-35 -15];
all_V_reset = [-70 -70];
    
all_E_L = all_V_th - cell_features_priors.dE_L-[0 6];
k = 0.07; % The stimulated intensity 
    
run('./LIF_simulation.m')



%% Visualization

figure(1)
plot(t_vect,voltage_traces{1,1} ,'col',[1,0,0,0.7],'Linewidth',2);
hold on;
ylim([-90,-10]);
xlim([0,max(t_vect)]);
plot(t_vect,voltage_traces{1,2} ,'col',[0,1,0,0.7],'Linewidth',2);
xlabel('Time (1/20 ms)');
ylabel('Presynaptic Membrane Potential');
hold off;


figure(2)
plot(presynaptic_events{1,1},presynaptic_amplitudes{1,1},'.','markers',20,'col', [1,0,0,0.8]);
hold on;
ylim([2,4]);
xlim([0,max(t_vect)]);
plot(presynaptic_events{1,2},presynaptic_amplitudes{1,2},'.','markers',20,'col',[0,1,0,0.8]);
xlabel('Time (1/20 ms)');
ylabel('Presynaptic events');
hold off;


figure(3)
plot(mpp_new(1).event_times,mpp_new(1).amplitudes,'.','markers',20,'col', [1,0,0,0.1]);
hold on;
ylim([2,4]);
xlim([0,max(t_vect)]);
xlabel('Time (1/20 ms)');
ylabel('Postsynaptic events');
hold off;


outflnm = strcat('../../Figures/Interpolation/FullModel_partI_');

%saveas(1,strcat(outflnm,'PrePotential','.jpg'));

%saveas(2,strcat(outflnm,'PreEvents','.jpg'));

%saveas(3,strcat(outflnm,'PostEvents','.jpg'));

%%
% Estimate the average intensity from the estimated parameters (by drawing
% multiple realizations)
run('./Expected_intensity.m');
%%
figure(4)
%plot(t_vect,KS_intensity{1,1},'col',[1,0,0,0.7],'Linewidth',2);
hold on;
ylim([-0.1,1]);
xlim([0,max(t_vect)]);
plot(t_vect,reshape(E_intensity{1,1},[length(t_vect) 1]),'--','col', [1,0,0,0.7]...
    ,'Linewidth',2);
hold on;
plot(t_vect(presynaptic_events{1,1}),-0.1,'.','markers',20,'col',[1,0,0,0.8]);
plot(t_vect,M_intensity,'-.','col', [1,0,0,0.7]...
    ,'Linewidth',2);

xlabel('Time (1/20 ms)');
ylabel('Monte-Carlo Intensity');
hold off;


figure(5)
plot(t_vect,KS_intensity{2,1},'col',[0,1,0,0.7],'Linewidth',2);
hold on;
plot(t_vect,reshape(E_intensity{2,1},[length(t_vect) 1]),'--','col', [0,1,0,0.7]...
    ,'Linewidth',2);

ylim([-0.1,1]);
xlim([0,max(t_vect)]);

plot(t_vect(presynaptic_events{1,2}),-0.1,'.','markers',20,'col',[0,1,0,0.8]);


xlabel('Time (1/20 ms)');
ylabel('Monte-Carlo Intensity');
hold off;

%'.','markers',20,

%saveas(4,strcat(outflnm,'Intensity_I','.jpg'));

%saveas(5,strcat(outflnm,'Intensity_II','.jpg'));

%% Run the forward-backward sampling algorithm
% Sampling using true presynaptic events 
rng(11142,'twister');

i_cell = 1;
assigned_events = presynaptic_events{1,i_cell};
run('./Gibbs_sampler_synaptic_failure.m');
V_seq_all = V_seq;
assigned_events_all = assigned_events;

% Create a random assignments:
i_cell = 1;
assigned_events = randsample(presynaptic_events{1,i_cell},3);

run('./Gibbs_sampler_synaptic_failure.m');
V_seq_partial = V_seq;
assigned_events_partial = assigned_events;
% Now rom a mixture
i_cell = 1;
assigned_events = randsample(mpp_new(1).event_times,3);

run('./Gibbs_sampler_synaptic_failure.m');
V_seq_mix = V_seq;
assigned_events_mix = assigned_events;


%% Visualization 

figure(6)
plot(t_vect(1:n_grid_time-1),V_seq_all(2:end,1),'col',[1,0,0,0.1],'Linewidth',2);
hold on;
ylim([-90,-10]);
xlim([0,max(t_vect)]);
for i = 1:n_sim_fb
    plot(t_vect(1:n_grid_time-1),V_seq_all(2:end,i),'--','col',[1,0,0,0.3], 'Linewidth',2);
end

plot(assigned_events_all,-80,'.','markers',20,'col', [1,0,0,0.1]);
hold on;
plot(t_vect,voltage_traces{1,1} ,'-','col',[0,0,1,0.8], 'Linewidth',2);
xlabel('Time (1/20 ms)');
ylabel('Sampled Voltage (true events)');
hold off;

figure(7)
plot(t_vect(1:n_grid_time-1),V_seq_partial(2:end,1),'col',[1,0,0,0.1],'Linewidth',2);
hold on;
ylim([-90,-10]);
xlim([0,max(t_vect)]);
for i = 1:n_sim_fb
    plot(t_vect(1:n_grid_time-1),V_seq_partial(2:end,i),'--','col',[1,0,0,0.3], 'Linewidth',2);
end
plot(t_vect,voltage_traces{1,1} ,'-','col',[0,0,1,0.8], 'Linewidth',2);
hold on;
plot(assigned_events_partial,-80,'.','markers',20,'col', [1,0,0,0.1]);

xlabel('Time (1/20 ms)');
ylabel('Sampled Voltage (w/ synaptic failure)');
hold off;


figure(8)
plot(t_vect(1:n_grid_time-1),V_seq_mix(2:end,1),'col',[1,0,0,0.1],'Linewidth',2);
hold on;
ylim([-90,-10]);
xlim([0,max(t_vect)]);
for i = 1:n_sim_fb
    plot(t_vect(1:n_grid_time-1),V_seq_mix(2:end,i),'--','col',[1,0,0,0.3], 'Linewidth',2);
end
plot(assigned_events_mix,-80,'.','markers',20,'col', [1,0,0,0.1]);
hold on;
plot(t_vect,voltage_traces{1,1} ,'-','col',[0,0,1,0.8], 'Linewidth',2);
xlabel('Time (1/20 ms)');
ylabel('Sampled Voltage (mixed events)');
hold off;

%saveas(6,strcat(outflnm,'Voltage_True','.jpg'));

%saveas(7,strcat(outflnm,'Voltage_Par','.jpg'));

%saveas(8,strcat(outflnm,'Voltage_Mix','.jpg'));

%% Draw a realization of the voltage from the probability 
% Assign the labels using simple probability ratio

classification = reshape(E_intensity{2,1} > E_intensity{1,1}, [length(t_vect) 1]);
figure(9)
plot(presynaptic_events{1,1},presynaptic_amplitudes{1,1},'.','markers',20,'col', [1,0,0,0.8]);
hold on;
ylim([2,4]);
xlim([0,max(t_vect)]);
plot(presynaptic_events{1,2},presynaptic_amplitudes{1,2},'.','markers',20,'col',[0,1,0,0.8]);

plot(t_vect+1,2.8*classification,'.','markers',20,'col',[0,1,0,0.8]);

plot(t_vect+1,2.8*(1-classification),'.','markers',20,'col',[1,0,0,0.8]);



xlabel('Time (1/20 ms)');
ylabel('Presynaptic events');
hold off;


%saveas(8,strcat(outflnm,'Intensity_III','.jpg'));
