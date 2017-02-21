%%
n_gibbs_sample = 100;
n_burnin = 500;
n_skip = 10;
% Initialize experiment conditions
n_cell = 20;
n_trial = 400;
n_sim_evoked =5; % evoked cells in each trial 
k_basic = 0.1;

exact_crossing = 0;

%% Run simulations 

rng(11111,'twister');

strength_seq = [0.1 0.1 2 2];
Vthre_std_seq = [0.1 4 0.1 4];
Vthreset_std_seq = [0.1 4 0.1 4];
tic
%for casenumber = 1:4
casenumber=3;
    connection_strength_stddev = strength_seq(casenumber); % size distribution
    Vthre_std = Vthre_std_seq(casenumber);% firing threshold distribution
    Vreset_std =Vthreset_std_seq(casenumber); % reset threshold distribution
    
    run('./LIF_simulation.m');
    run('./Simulation_crude.m');
    run('./Expected_intensity.m');
    
%     tstart=toc;
%     n_trial_update = 80; % number of trials to update in each iteration
%     run('./Simulation_medium.m');
%     tend=toc;
%     t_delta = tend-tstart;
%     flnm=strcat('../../Data/Full_sim/Case', num2str(casenumber),'_E_joint_minibatch.mat');
%     save(flnm,'t_delta','mpp_new','sigma_samples','gamma_samples','mu_samples','all_connected', 'all_amplitudes',...
%         'assignments_samples');
%     
%      sigma_unknown=1;
%     n_trial_update = 400; % number of trials to update in each iteration
%     tstart=toc;
%     run('./Simulation_medium.m');
%      tend=toc;
%       t_delta = tend-tstart;
%    
%     flnm=strcat('../../Data/Full_sim/Case', num2str(casenumber),'_E_joint.mat');
%     save(flnm,'t_delta','mpp_new','sigma_samples','gamma_samples','mu_samples','all_connected', 'all_amplitudes',...
%         'assignments_samples');

    sigma_unknown=0;

%     n_trial_update = 400; % number of trials to update in each iteration    
%     tstart=toc;
%     run('./Simulation_medium.m');
%      tend=toc;
%       t_delta = tend-tstart;
%    
%     flnm=strcat('../../Data/Full_sim/Case', num2str(casenumber),'_E_joint_known.mat');
%     save(flnm,'t_delta','mpp_new','sigma_samples','gamma_samples','mu_samples','all_connected', 'all_amplitudes',...
%         'assignments_samples');
    
         n_trial_update = 80;
          sigma_unknown=0;
    tstart=toc;
    run('./Simulation_medium.m');
     tend=toc;
      t_delta = tend-tstart;
   
    flnm=strcat('../../Data/Full_sim/Case', num2str(casenumber),'_E_joint_known_minibatch.mat');
    save(flnm,'t_delta','mpp_new','sigma_samples','gamma_samples','mu_samples','all_connected', 'all_amplitudes',...
        'assignments_samples');
    
         n_trial_update = 80;
          sigma_unknown=0;
    tstart=toc;
    run('./Simulation_integral.m');
     tend=toc;
      t_delta = tend-tstart;
   
    flnm=strcat('../../Data/Full_sim/Case', num2str(casenumber),'_E_joint_known_minibatch_int.mat');
    save(flnm,'t_delta','mpp_new','sigma_samples','gamma_samples','mu_samples','all_connected', 'all_amplitudes',...
        'assignments_samples');
    
    
%end
%% Evaluate the posterior

%for casenumber=1:4
casenumber = 3;

% flnm=strcat('../../Data/Full_sim/Case', num2str(casenumber),'_E_joint.mat');
% load(flnm)
% mean_gamma = mean(gamma_samples,1);
% mean_amp = mean(mu_samples,1);
% mean_sigma= mean(sigma_samples,1);
% comptime=t_delta;
% 
% flnm=strcat('../../Data/Full_sim/Case', num2str(casenumber),'_E_joint_known.mat');
% load(flnm)
% mean_gamma_known = mean(gamma_samples,1);
% mean_amp_known = mean(mu_samples,1);
% mean_sigma_known= mean(sigma_samples,1);
% comptime_known=t_delta;

% flnm=strcat('../../Data/Full_sim/Case', num2str(casenumber),'_E_joint_minibatch.mat');
% load(flnm)
% mean_gamma_mini = mean(gamma_samples,1);
% mean_amp_mini = mean(mu_samples,1);
% mean_sigma_mini= mean(sigma_samples,1);
% t_delta


flnm=strcat('../../Data/Full_sim/Case', num2str(casenumber),'_E_joint_known_minibatch.mat');
load(flnm)
mean_gamma_known_mini = mean(gamma_samples,1);
mean_amp_known_mini = mean(mu_samples,1);
mean_sigma_known_mini= mean(sigma_samples,1);
comptime_known_mini=t_delta;

flnm=strcat('../../Data/Full_sim/Case', num2str(casenumber),'_E_joint_known_minibatch_int.mat');
load(flnm)
mean_gamma_known_mini_int = mean(gamma_samples,1);
mean_amp_known_mini_int = mean(mu_samples,1);
mean_sigma_known_mini_int= mean(sigma_samples,1);
comptime_known_mini_int=t_delta;
mean_gamma_known_mini_int

flnm=strcat('../../Data/Full_sim/Case', num2str(casenumber),'Crude.mat');
load(flnm)

n_connected = sum(all_connected);
n_disconnected = n_cell - n_connected;

figure(1)
% plot(1:n_connected,mean_gamma(find(all_connected)),'.','markers',20,'col',[1,0,0,0.8]);
% hold on;
% plot(n_connected+(1:n_disconnected),mean_gamma(find(all_connected==0)),'.','markers',10,'col',[1,0,0,0.8]);
% hold on;
% 
% % With sigma assumed to be known 
% plot(1:n_connected,mean_gamma_known(find(all_connected)),'.','markers',20,'col',[0,0,1,0.8]);
% hold on;
% plot(n_connected+(1:n_disconnected),mean_gamma_known(find(all_connected==0)),'.','markers',10,'col',[0,0,1,0.8]);
% hold on;

% With sigma assumed to be known and minibatch
plot(1:n_connected,mean_gamma_known_mini(find(all_connected)),'.','markers',20,'col',[0,1,0,0.8]);
hold on;
plot(n_connected+(1:n_disconnected),mean_gamma_known_mini(find(all_connected==0)),'.','markers',10,'col',[0,1,0,0.8]);
hold on;

% With sigma assumed to be known and minibatch
plot(1:n_connected,mean_gamma_known_mini_int(find(all_connected)),'.','markers',20,'col',[0,0,1,0.8]);
hold on;
plot(n_connected+(1:n_disconnected),mean_gamma_known_mini_int(find(all_connected==0)),'.','markers',10,'col',[0,0,1,0.8]);
hold on;

% The crude estimates 
plot(1:n_connected,overall_connectivity(find(all_connected)),'o','markers',10,'col',[1,0,0,0.8]);
hold on;
plot(n_connected+(1:n_disconnected),overall_connectivity(find(all_connected==0)),'o','markers',10,'col',[1,0,0,0.8]);
hold on;

%bar(-5,0.5*comptime/comptime,'r')

%bar(-3,0.5*comptime_known/comptime,'b')

%bar(-1,0.5*comptime_known_mini/comptime,'g')

xticks([-5 -3 -1 1 5 10 15 20])
xticklabels({'A','B','C', '1', '5', '10', '15', '20'})


line([0 0],[1 0]); 
line([n_connected+0.5 n_connected+0.5],[8 0]); 
ylim([0,1]);
xlim([-6,20]);
xlabel('Cells');
ylabel('Posterior mean of Gamma');
hold off;

% 
figure(2)
% plot(1:n_connected,mean_amp(find(all_connected)),'.','markers',20,'col',[1,0,0,0.8]);
% hold on;
 plot(1:n_connected,all_amplitudes(all_connected),'*','markers',10,'col',[0,0,0,0.8]);
 hold on;
plot(1:n_connected,overall_mark(all_connected),'o','markers',10,'col',[0,0,0,0.8]);
hold on;
plot(1:n_connected,mean_amp_known_mini(find(all_connected)),'.','markers',20,'col',[0,1,0,0.8]);
hold on;
plot(1:n_connected,mean_amp_known_mini_int(find(all_connected)),'.','markers',20,'col',[0,0,1,0.8]);

% plot(1:n_connected,mean_amp_known(find(all_connected)),'.','markers',20,'col',[0,0,1,0.8]);
% hold on;
% plot(n_connected+(1:n_disconnected),mean_amp(find(all_connected==0)),'.','markers',20,'col',[1,0,0,0.8]);
% hold on;
 plot(n_connected+(1:n_disconnected),all_amplitudes(all_connected==0),'*','markers',10,'col',[0,0,0,0.8]);
 hold on;
plot(n_connected+(1:n_disconnected),overall_mark(all_connected==0),'o','markers',10,'col',[0,0,0,0.8]);
hold on;
plot(n_connected+(1:n_disconnected),mean_amp_known_mini(find(all_connected==0)),'.','markers',20,'col',[0,1,0,0.8]);
hold on;
% plot(n_connected+(1:n_disconnected),mean_amp_known(find(all_connected==0)),'.','markers',20,'col',[0,0,1,0.8]);
plot(n_connected+(1:n_disconnected),mean_amp_known_mini_int(find(all_connected==0)),'.','markers',20,'col',[0,0,1,0.8]);

line([n_connected+0.5 n_connected+0.5],[8 0]); 
ylim([0,8]);
xlim([0,20]);
xlabel('Cells');
ylabel('Posterior mean of mu');
hold off;


%saveas(1,strcat('../../Figures/Full_model/S',num2str(casenumber),'Gamma.jpg'))
%saveas(2,strcat('../../Figures/Full_model/S',num2str(casenumber),'mu.jpg'))
%end
%% Draw the assignments:
casenumber=1;
flnm=strcat('../../Data/Full_sim/Case', num2str(casenumber),'_E_joint_known_minibatch.mat');
load(flnm)

cell_id = 10;
assignments_collapse = [mpp_new.assignments];
true_events_this_cell = find(assignments_collapse==cell_id);
amplitudes_collapse = [mpp_new.amplitudes];

true_events_amp = amplitudes_collapse(true_events_this_cell);
gibbs_assignments_one = [assignments_samples{70}{:}];
true_events_gibbs = gibbs_assignments_one(true_events_this_cell);
list_cells = unique(true_events_gibbs);


all_events_amp = cell(length(list_cells),1);
counts=cell(length(list_cells),1);
binrng=(1:10)*(max(true_events_amp) - min(true_events_amp))/9 +min(true_events_amp);

for i_cell = 1:length(list_cells)
    all_events_amp{i_cell} = true_events_amp(find(true_events_gibbs==list_cells(i_cell)));
    if i_cell == 1
        counts{i_cell} = histc(all_events_amp{i_cell}, binrng);     
    else
        counts{i_cell} = counts{i_cell-1}+histc(all_events_amp{i_cell}, binrng);     
    end
end

colorset =jet(length(list_cells));
%colorset = [colorset ones(length(list_cells),1)*0.7];
figure(3)
for i_cell = 1:length(list_cells)
    bar(binrng, counts{1+length(list_cells)-i_cell},'FaceColor',...
        colorset(1+length(list_cells)-i_cell,:));
    hold on
end
legend(cellstr(num2str(fliplr(list_cells)', '%-d')));

xlabel('Amplitudes');
ylabel('Counts');
title(strcat('Assignments of events from Cell  ', num2str(cell_id),' from one Gibbs sample'));
saveas(3,strcat('../../Figures/Full_model/Assignments_1.jpg'))

%%

casenumber=1;
flnm=strcat('../../Data/Full_sim/Case', num2str(casenumber),'_E_joint_known_minibatch.mat');
load(flnm)

assignments_collapse = [mpp_new.assignments];
events_index = [89 219];

true_index = assignments_collapse(events_index);
assigned_values= zeros([100 2]);
for i_gibbs = 1:100
    gibbs_assignments_one = [assignments_samples{i_gibbs}{:}];
    assigned_values(i_gibbs,:)=gibbs_assignments_one(events_index);
end
true_index

figure(4)
plot(1:100,assigned_values(:,1),'.','markers',10,'col',[0,0,1,0.8]);
hold on;
plot(1:100,assigned_values(:,1),'col',[0,0,1,0.5]);
hold on;

ylim([0,20]);
xlim([0, 100]);
xlabel('Index of samples');
ylabel('Assignments');

title(strcat('Assignments of events from Cell  ', num2str(true_index(1)),''));
saveas(4,strcat('../../Figures/Full_model/Assignments_gibbs_1.jpg'))


figure(5)
plot(1:100,assigned_values(:,2),'.','markers',10,'col',[0,0,1,0.8]);
hold on;
plot(1:100,assigned_values(:,2),'col',[0,0,1,0.5]);
hold on;

ylim([0,20]);
xlim([0, 100]);
xlabel('Index of samples');
ylabel('Assignments');

title(strcat('Assignments of events from Cell  ', num2str(true_index(2)),''));
saveas(5,strcat('../../Figures/Full_model/Assignments_gibbs_2.jpg'))
