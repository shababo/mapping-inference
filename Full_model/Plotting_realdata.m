%% Plot the cells with estimated gammas 
%% Load the gamma estimates and data 
load('./Environments/6_3_s2c3_mpp_and_stim_data.mat')
load('final_fits_6_3_single_reduced.mat');%'gamma_all','gain_all');
load('initial_fits_6_3_single_reduced.mat'); %'gamma_initial');
outflnm =strcat('./Figures/Plots_6_3_single_reduced');


% load('./Environments/6_5_s2c1_mpp_and_stim_data.mat')
% load('final_fits_6_5_single.mat');%'gamma_all','gain_all');
% load('initial_fits_6_5_single.mat'); %'gamma_initial');
% outflnm =strcat('./Figures/Plots_6_5_single');



%%
local_locations= cell_locs;
Z_dense =target_locs;
powers_trials = stim_pow;
locations_trials=target_inds;

%% Initial fits
figure(1)


for i = 1:n_cell_local
    if gamma_initial_all(i)>0.05
temp2 = scatter(local_locations(i,2)+151,local_locations(i,1)+151,...
    200*gamma_initial_all(i),'MarkerEdgeColor','r','MarkerFaceColor','r',...
    'MarkerFaceAlpha',(local_locations(i,3)/max(local_locations(:,3))));
    end
hold on;
end
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
%%

figure(2)
temp2 = scatter(local_locations(gain_initial_all>0.001,2)+151,local_locations(gain_initial_all>0.001,1)+151,...
    4000*gain_initial_all(gain_initial_all>0.001));
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


%% First fit 
figure(3)

for i = 1:n_cell_local
    if gamma_first_fit(i)>0.05
temp2 = scatter(local_locations(i,2)+151,-local_locations(i,1)-151,...
    200*gamma_first_fit(i),'MarkerEdgeColor','r','MarkerFaceColor','r',...
    'MarkerFaceAlpha',(local_locations(i,3)/max(local_locations(:,3))));
    end
hold on;
end


% Locate  one isolated cell
% this_cell = find(local_locations(:,2) < 130 & local_locations(:,2) > 100 & local_locations(:,1) < 50 & local_locations(:,1) > 30);
% gamma_first_fit(this_cell)
% this_cell=131;
% temp2 = scatter(local_locations(this_cell,2)+151,local_locations(this_cell,1)+151,...
%     200*gamma_first_fit(this_cell),'MarkerEdgeColor','g','MarkerFaceColor','g',...
%     'MarkerFaceAlpha',0.7);


temp1 = scatter(local_locations(:,2)+151,-local_locations(:,1)-151,...
    2);
set(temp1,'MarkerEdgeColor','g','MarkerFaceColor','g');
alpha(temp1,1);
hold on;


set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])

xlim([-20,313]);
ylim([-313,20]);
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
 saveas(3,strcat(outflnm,'Gamma_first','.jpg'));


%% Final fits 

figure(4)

for i = 1:n_cell_local
    if gamma_final_all(i)>0.05
temp2 = scatter(local_locations(i,2)+151,local_locations(i,1)+151,...
    200*gamma_final_all(i),'MarkerEdgeColor','r','MarkerFaceColor','r',...
    'MarkerFaceAlpha',(local_locations(i,3)/max(local_locations(:,3))));
    end
hold on;
end

temp1 = scatter(local_locations(:,2)+151,local_locations(:,1)+151,...
    2);
set(temp1,'MarkerEdgeColor','g','MarkerFaceColor','g');
alpha(temp1,1);
hold on;


set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])

xlim([-20,313]);
ylim([-20,313]);
% axis off;
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

% saveas(4,strcat(outflnm,'Gamma_final','.jpg'));





