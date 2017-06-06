%
addpath(genpath('c:/Users/Shizhe/Documents/GitHub/psc-detection'),...
    genpath('C:/Users/Shizhe/Documents/GitHub/mapping-inference'),...
    genpath('C:/Users/Shizhe/Documents/GitHub/mapping-core'));
%% Figure 11: single spot cell maps
current_wd = pwd;
%%
% current_wd = 'C:/Users/Shizhe/Dropbox/Presentations/2017-5-30 SAND8/Code';
%%
run('./Experiments_I.m');

%%

cd(current_wd);

    gamma_initial =zeros(n_cell_local,1);
    gamma_initial(stimulated_cells)=gamma_path(:,end);
save('initial_fits.mat','gamma_initial');

%% Summary stat:
% sum([mpp.assignments]==0) %1560
%     300*n_trial*background_rate %1552
%     
%     
%     length([mpp.assignments]) %2342
%  length([mpp_copy(1:n_trial).times]) %2269
%  

%% Draw the cell shape map:

shape_map = zeros(301,301);

for i_cell = 1:n_cell_local
    cell_center = local_locations(i_cell,1:2);
    cell_center=round(cell_center)+[151 151];
    for i= 11:71
        i_map=cell_center(1) + i-41;
        
        
        if i_map> 0 & i_map<302
            for j = 11:71
                j_map=cell_center(2)+j-41;
                if j_map>0 & j_map<302
                    if isnan(l23_cells_for_sim(local_shape_gain(i_cell)).shape(i,j,91))
                        l23_cells_for_sim(local_shape_gain(i_cell)).shape(i,j,91)=0;
                    end
                    shape_map(i_map,j_map)=...
                        max(shape_map(i_map,j_map), l23_cells_for_sim(local_shape_gain(i_cell)).shape(i,j,91));
                end
            end
        end
    end
end
% 
% %%
% cell_loc = zeros(0,2);
% for i_map = 1:301
%     for j_map = 1:301
%          if shape_map(i_map,j_map)>0.6
%             cell_loc = [cell_loc; [i_map j_map]]; 
%          end
%     end
% end
%%
figure(1)
heatmap(shape_map','Colormap','redgreencmap');
hold on;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])

temp = scatter(local_locations(:,1)+151,local_locations(:,2)+151,...
    15,'fill','k');
hold on;
set(temp,'MarkerFaceColor','k');
alpha(temp,0.5);

outflnm =strcat('../Figures/Results/');
saveas(1,strcat(outflnm,'CellShape','.jpg'));
%%

figure(2)

temp = scatter(local_locations(local_gamma>0.05,2)+151,local_locations(local_gamma>0.05,1)+151,...
    300*local_gamma(local_gamma>0.05));
hold on;
set(temp,'MarkerFaceColor','b');
alpha(temp,0.5);


temp2 = scatter(local_locations(gamma_initial>0.05,2)+151,local_locations(gamma_initial>0.05,1)+151,...
    300*gamma_ini(gamma_ini>0.05));
set(temp2,'MarkerFaceColor','r');
alpha(temp2,0.5);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])

xlim([0,301]);
ylim([0,301]);
axis off;
outflnm =strcat('../Figures/Results/');
saveas(2,strcat(outflnm,'Gamma_ini','.jpg'));



%%
run('./Experiments_II.m');
%%
% cd(current_wd);

gamma_current = gamma_fits(:,num_iter-1);
gamma_all = zeros(n_cell_local,1);
gamma_all(stimulated_cells)=gamma_current;
local_gain = [l23_cells_for_sim(local_shape_gain).optical_gain]';
gain_all = zeros(n_cell_local,1);
gain_all(stimulated_cells)=gain_current;



save('final_fits_Set2.mat','local_gamma','local_gain','gamma_all','gain_all');

%%
post_shape_map = zeros(301,301);

for i_cell = 1:n_cell_local
    cell_center = local_locations(i_cell,1:2);
    cell_center=round(cell_center)+[151 151];
    if gamma_all(i_cell)>0.3
    for i= 11:71
        i_map=cell_center(1) + i-41;
        
        
        if i_map> 0 & i_map<302
            for j = 11:71
                j_map=cell_center(2)+j-41;
                if j_map>0 & j_map<302
                    if isnan(l23_cells_for_sim(local_shape_gain(i_cell)).shape(i,j,91))
                        l23_cells_for_sim(local_shape_gain(i_cell)).shape(i,j,91)=0;
                    end
                    post_shape_map(i_map,j_map)=...
                        max(post_shape_map(i_map,j_map), l23_cells_for_sim(local_shape_gain(i_cell)).shape(i,j,91));
                end
            end
        end
    end
    end
end

% 
% %%
% cell_loc = zeros(0,2);
% for i_map = 1:301
%     for j_map = 1:301
%          if shape_map(i_map,j_map)>0.6
%             cell_loc = [cell_loc; [i_map j_map]]; 
%          end
%     end
% end
%%
cd(current_wd);
figure(3)
heatmap(post_shape_map','Colormap','redgreencmap');
hold on;

temp = scatter(local_locations(gamma_all>0.1,1)+151,local_locations(gamma_all>0.1,2)+151,...
    30,'fill','r');
hold on;
set(temp,'MarkerFaceColor','r');
alpha(temp,0.5);


set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])

outflnm =strcat('../Figures/Results/');
saveas(3,strcat(outflnm,'CellShape_Post','.jpg'));

%%
load('final_fits.mat')

figure(4)
temp = scatter(local_locations(local_gamma>0.05,2)+151,local_locations(local_gamma>0.05,1)+151,...
    200*local_gamma(local_gamma>0.05));
hold on;
set(temp,'MarkerEdgeColor','b','MarkerFaceColor','b');
alpha(temp,0.5);


temp1 = scatter(local_locations(:,2)+151,local_locations(:,1)+151,...
    2);
set(temp1,'MarkerEdgeColor','g','MarkerFaceColor','g');
alpha(temp1,1);
hold on;


temp2 = scatter(local_locations(gamma_all>0.05,2)+151,local_locations(gamma_all>0.05,1)+151,...
    200*gamma_all(gamma_all>0.05));
set(temp2,'MarkerEdgeColor','r','MarkerFaceColor','r');
alpha(temp2,0.3);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])

xlim([-20,313]);
ylim([-20,313]);
axis off;
% Add legends to the plot 
temp5 = scatter([-9 -9 -9],[300 280 260],...
    190*[1 0.5 0.2]);
set(temp5,'MarkerEdgeColor','r','MarkerFaceColor','r');
alpha(temp5,0.3);



txt4 = '\gamma = 1';
txt5 = '\gamma = 0.5';
txt6 = '\gamma = 0.2';

text(4,300,txt4)

text(4,280,txt5)

text(4,260,txt6)



temp4 = scatter([65 65],[300 260],...
    200*[1 1]);
set(temp4,'MarkerEdgeColor','b','MarkerFaceColor','b');
alpha(temp4,0.3);

temp3 = scatter([65 65],[280 260],...
    190*[1 1]);
set(temp3,'MarkerEdgeColor','r','MarkerFaceColor','r');
alpha(temp3,0.3);


txt1 = 'Truth';
txt2 = 'Estimate';
txt3 = 'Good fit';

text(78,300,txt1)

text(78,280,txt2)

text(78,260,txt3)


rectangle('Position',[-20 246 155 67],'LineStyle','--')


outflnm =strcat('../Figures/Results/');
saveas(4,strcat(outflnm,'Gamma_final','.jpg'));

%%
%Draw the plot for real data
%load('final_fits_data_new.mat')

figure(5)
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


outflnm =strcat('../Figures/Results/');
saveas(5,strcat(outflnm,'Gamma_final_data_new','.jpg'));

%% Check the amount of stimuli received by the connected cells
high_gamma =sum( stimuli_size_local(:,gamma_all==1));
low_gamma =sum( stimuli_size_local(:,gamma_all<1 & gamma_all>0.2));
no_gamma =sum( stimuli_size_local(:,gamma_all<0.2));

figure(1)
plot(1:length(high_gamma),high_gamma,'.','MarkerSize',20,'col','r')
hold on;
plot(length(high_gamma)+(1:length(low_gamma)),low_gamma,'.','MarkerSize',20,'col','b')
hold off;
ylabel('Total stimulations')
%%
stim_gamma =sum( stimuli_size_local(:,:));
figure(1)
plot(gamma_all,stim_gamma,'.','MarkerSize',20,'col','r')
hold on;
ylabel('Total stimulations')
xlabel('Gamma fits')




