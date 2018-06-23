%% Simulate a two-neuron system with delays:
addpath(genpath('../../../GitHub/mapping-inference'));
bg_rate=1e-4;
%% Use real data (z-axis and xy )
z_path='../Environments/z_data.mat';
xy_path='../Environments/new_xy_data.mat';

axis_list={'x' 'y' 'z'};

for i_ax = 1:3 % choose the axis
arg1=1; % random seed, set it to others for random initialization
% i_ax=2;
%% Scale the current (devide by power and then by the max curr)

%% Reformat data into neurons
target_axis=axis_list{i_ax};
[neurons] = map_result_to_neuron(target_axis, z_path, xy_path);

%% Preprocessing:
% 1) find modes
% 2) estimate mean function
% 3) estimate variance function

params=struct;
if i_ax==1
    params.boundary= 35; % two-side
    params.buffer=10;
    params.mag=1;
    params.tau=49;
elseif i_ax==2
    params.boundary= 25; % two-side
    params.buffer=10;
    params.mag=1;
    params.tau=49;
else
    params.boundary= 60; % two-side
    params.buffer=10;
    params.mag=1;
    params.tau=400;
end
params.symmetric=false;
[mean_params, var_params,neurons]=pre_processing(neurons,params);

% Set initial values:
rng(arg1)
[variational_params, prior_params,inference_params]=initialization_for_shape(i_ax,...
    neurons,mean_params,var_params);
% Use VI to learn the tau for GP kernel

[parameter_hist,lklh_hist,elbo_hist]=fit_shape_VI(neurons,variational_params,prior_params,...
    inference_params);
%
file_path = ['./Data/New_Real_Data_' target_axis 'ini' num2str(arg1) '.mat'];
save(file_path,'neurons', 'parameter_hist','lklh_hist','elbo_hist',...
    'var_params','mean_params')



%% Plotting:
% fig_name=['./Figures/OLplot_' target_axis];
% 
% i_cell = 5;
% 
% fig=figure(1);
% subplot(1,2,1)  
% 
% scatter(neurons(i_cell).stim_grid,neurons(i_cell).raw_current,...
%     'MarkerFaceColor','r')
% hold on;
% % plot(neurons(i_cell).unique_grid,neurons(i_cell).avg_current,'Color','r')
% % hold on;
% 
% title([target_axis '-axis, cell ' num2str(i_cell)])
% xlabel('Stim location')
% ylabel('Raw current')
%  
% 
% subplot(1,2,2) 
% scatter(neurons(i_cell).stim_grid,neurons(i_cell).scaled_current,...
%     'MarkerFaceColor','k')
% hold on;
% 
% scatter(neurons(i_cell).stim_grid,neurons(i_cell).power/(max(neurons(i_cell).power)),...
%     'MarkerFaceColor','g')
% 
%  xlabel('Stim location')
%  ylabel('Scaled current')
%  title([target_axis '-axis, cell ' num2str(i_cell)])
% 
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 10 4];
% saveas(fig,[fig_name  'scale.png'])
% close(fig)
% %% Visualize the scaled variance:
% 
% fig=figure(2);
% subplot(1,2,1)
% 
% scatter(neurons(i_cell).unique_grid,neurons(i_cell).raw_sq_deviation,...
%     'MarkerFaceColor','r')
% hold on;
% 
% title([target_axis '-axis, cell ' num2str(i_cell)])
% xlabel('Stim location')
% ylabel('Squared deviation')
% 
% 
% subplot(1,2,2)
% scatter(neurons(i_cell).unique_grid,(neurons(i_cell).sq_deviation)/(neurons(i_cell).scale),...
%     'MarkerFaceColor','k')
% hold on;
% 
% 
% xlabel('Stim location')
% ylabel('Scaled sq. dev')
% title([target_axis '-axis, cell ' num2str(i_cell)])
% 
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 10 4];
% saveas(fig,[fig_name  'scale_variance.png'])
% close(fig)
% 
% 
% %% Visualize the mode finding performance:
% 
% n_cell = length(neurons);
% colors=lines(n_cell);
% 
% fig=figure(3);
% subplot(1,2,1)
% for i_cell = 1:n_cell
% scatter(neurons(i_cell).stim_grid,neurons(i_cell).scaled_current,...
%     'MarkerFaceColor',colors(i_cell,:))
% hold on;
% end
% title([target_axis '-axis'])
% ylabel('Normalized current')
% 
% xlabel('Stim location')
% 
% subplot(1,2,2)
% for i_cell = 1:n_cell
% scatter(neurons(i_cell).adjusted_grid,neurons(i_cell).scaled_current,...
%     'MarkerFaceColor',colors(i_cell,:))
% hold on;
% end
% hold on;
% 
% 
% xlabel('Adjusted stim location')
% ylabel('Normalized current')
% title([target_axis '-axis; mean shift: ' num2str(round(mean([neurons(:).initial_shift] ),1) )])
% 
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 10 4];
% saveas(fig,[fig_name  'mode.png'])
% close(fig)
% 
% %% Draw the predictive mean function:
% 
% fig=figure(4);
% 
% subplot(1,2,1)
% scatter(mean_params.data.X,mean_params.data.Y,...
%     'MarkerFaceColor','r')
% hold on;
% plot(mean_params.grid,mean_params.values,'linewidth',2,'Color','k')
% title([target_axis '-axis ' 'mean function'])
% xlabel('Stim location')
% 
% 
% subplot(1,2,2)
% plot(mean_params.grid,mean_params.prior_var,'linewidth',2)
% 
% title([target_axis '-axis prior variance'])
% ylim([min(mean_params.prior_var) max(mean_params.prior_var)+0.2 ])
% xlabel('Stim location')
% % ylabel('Normalized current')
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 10 4];
% saveas(fig,[fig_name  'mean.png'])
% close(fig)
% 
% 
% 
% %% Draw the predictive variance function:
% 
% fig=figure(5);
% 
% subplot(1,2,1)
% scatter(var_params.data.X,var_params.data.Y,...
%     'MarkerFaceColor','r')
% hold on;
% plot(var_params.grid,var_params.values,'linewidth',2,'Color','k')
% title([target_axis '-axis ' 'variance function'])
% xlabel('Stim location')
% % ylabel('Normalized current')
% 
% 
% subplot(1,2,2)
% plot(var_params.grid,var_params.prior_var,'linewidth',2)
% title([target_axis '-axis prior variance'])
% ylim([min(mean_params.prior_var) max(mean_params.prior_var)+0.2 ])
% xlabel('Stim location')
% % ylabel('Normalized current')
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 10 4];
% saveas(fig,[fig_name  'var.png'])
% close(fig)
end
