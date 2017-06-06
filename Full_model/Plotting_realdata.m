%% Plot the cells with estimated gammas 
%% Load the gamma estimates and data 
% load('./Environments/6_3_s2c3_mpp_and_stim_data.mat')
% load('final_fits_6_3.mat');%'gamma_all','gain_all');
% load('initial_fits_6_3.mat'); %'gamma_initial');
% outflnm =strcat('./Figures/Plots_6_3');


load('./Environments/6_5_s2c1_mpp_and_stim_data.mat')
load('final_fits_6_5.mat');%'gamma_all','gain_all');
load('initial_fits_6_5.mat'); %'gamma_initial');
outflnm =strcat('./Figures/Plots_6_5');



%%
local_locations= cell_locs;
Z_dense =target_locs;
powers_trials = stim_pow;
locations_trials=target_inds;

%% Initial fits
figure(1)
temp2 = scatter(local_locations(gamma_initial>0.05,2)+151,local_locations(gamma_initial>0.05,1)+151,...
    200*gamma_initial(gamma_initial>0.05));
set(temp2,'MarkerEdgeColor','r','MarkerFaceColor','r');
alpha(temp2,0.3);
hold on;

temp1 = scatter(local_locations(:,2)+151,local_locations(:,1)+151,...
    2);
set(temp1,'MarkerEdgeColor','g','MarkerFaceColor','g');
alpha(temp1,1);
hold on;


set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])

xlim([-20,313]);
ylim([-20,313]);
axis off;
% Add legends to the plot 
temp5 = scatter([0 0 0],[300 280 260],...
    200*[1 0.5 0.2]);
set(temp5,'MarkerEdgeColor','r','MarkerFaceColor','r');
alpha(temp5,0.3);


txt4 = '\gamma = 1';
txt5 = '\gamma = 0.5';
txt6 = '\gamma = 0.2';

text(12,300,txt4)

text(12,280,txt5)

text(12,260,txt6)

rectangle('Position',[-12 250 80 63],'LineStyle','--')

saveas(1,strcat(outflnm,'Gamma_initial','.jpg'));


%% Final fits 

figure(2)
temp2 = scatter(local_locations(gamma_all>0.05,2)+151,local_locations(gamma_all>0.05,1)+151,...
    200*gamma_all(gamma_all>0.05));
set(temp2,'MarkerEdgeColor','r','MarkerFaceColor','r');
alpha(temp2,0.3);
hold on;


temp1 = scatter(local_locations(:,2)+151,local_locations(:,1)+151,...
    2);
set(temp1,'MarkerEdgeColor','g','MarkerFaceColor','g');
alpha(temp1,1);
hold on;


set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])

xlim([-20,313]);
ylim([-20,313]);
axis off;
% Add legends to the plot 
temp5 = scatter([0 0 0],[300 280 260],...
    200*[1 0.5 0.2]);
set(temp5,'MarkerEdgeColor','r','MarkerFaceColor','r');
alpha(temp5,0.3);


txt4 = '\gamma = 1';
txt5 = '\gamma = 0.5';
txt6 = '\gamma = 0.2';

text(12,300,txt4)

text(12,280,txt5)

text(12,260,txt6)

rectangle('Position',[-12 250 80 63],'LineStyle','--')

saveas(2,strcat(outflnm,'Gamma_final','.jpg'));





