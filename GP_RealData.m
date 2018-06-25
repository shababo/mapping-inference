%% Simulate a two-neuron system with delays:
addpath(genpath('../../../GitHub/mapping-inference'));
bg_rate=1e-4;
%% Use real data (z-axis and xy )
z_path='../Environments/z_data.mat';
xy_path='../Environments/new_xy_data.mat';

axis_list={'x' 'y' 'z'};

for i_ax = 1:3 % choose the axis
% arg1=1; % random seed, set it to others for random initialization
i_ax=3;
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
    params.boundary= 45; % two-side
    params.buffer=10;
    params.mag=1;
    params.tau=49;
elseif i_ax==2
    params.boundary= 40; % two-side
    params.buffer=10;
    params.mag=1;
    params.tau=49;
else
    params.boundary= 100; % two-side
    params.buffer=10;
    params.mag=1;
    params.tau=400;
end
params.power_scaling=false;
params.symmetric=false;
params.sigma_current_model=true;

[mean_params, var_params,neurons]=pre_processing(neurons,params);

% Set initial values:
rng(arg1)
[variational_params, prior_params,inference_params]=initialization_for_shape(i_ax,...
    neurons,mean_params,var_params);
% Use VI to learn the tau for GP kernel

% [parameter_hist,lklh_hist,elbo_hist]=fit_shape_VI(neurons,variational_params,prior_params,...
%     inference_params);
% 
% file_path = ['./Data/New_Real_Data_' target_axis 'ini' num2str(arg1) '.mat'];
% save(file_path,'neurons', 'parameter_hist','lklh_hist','elbo_hist',...
%     'var_params','mean_params')


% output_plot=false;
% 
% % Plotting:
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
%  if output_plot
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 10 4];
% saveas(fig,[fig_name  'scale.png'])
% close(fig)
%  end
% %% Visualize the scaled variance:
% 
% fig=figure(2);
% 
% subplot(1,2,1)
% n_cell = length(neurons);
% colors=lines(n_cell);
% noise_var=[];mean_current =[];
%        
% for i_cell =1:n_cell
%     noise_var = [noise_var neurons(i_cell).raw_sq_deviation'];
%     mean_current = [mean_current neurons(i_cell).avg_current'];
%     
% end
% fitted_Y  = neurons(i_cell).slope* mean_current;
% 
% % fitted_sigma=[neurons(:).noise_sigma];
% % mean_raw_current = [neurons(:).mean_raw_current];
% for i_cell = 1:n_cell
%     scatter(neurons(i_cell).avg_current,sqrt(neurons(i_cell).raw_sq_deviation),'MarkerFaceColor',colors(i_cell,:))
%     hold on;
% end
% plot(reshape(mean_current,[size(mean_current,1)*size(mean_current,2) 1]),...
%     reshape(fitted_Y,[size(fitted_Y,1)*size(fitted_Y,2) 1]),...
%     'Color','k','LineWidth',3)
% hold on;
% title(['Current ' target_axis '-axis'])
% ylabel('Standard deviation')
% xlabel('Mean')
% 
% 
% 
% 
% 
% 
% subplot(1,2,2)
% scatter(neurons(i_cell).stim_grid,neurons(i_cell).scaled_current,...
%     'MarkerFaceColor','r')
% hold on;
% 
% scatter(neurons(i_cell).stim_grid,neurons(i_cell).noise_sigma,...
%     'MarkerFaceColor','k')
% 
% title([target_axis '-axis, cell ' num2str(i_cell)])
% xlabel('Stim location')
% ylabel('Normalized current and s.d.')
% 
% 
%  if output_plot
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 10 4];
% saveas(fig,[fig_name  'scale_variance.png'])
% close(fig)
% 
%  end
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
% 
%  if output_plot
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 10 4];
% saveas(fig,[fig_name  'mode.png'])
% close(fig)
%  end
% %% Draw the predictive mean function:
% 
% fig=figure(4);
% 
% scatter(mean_params.data.X,mean_params.data.Y,...
%     'MarkerFaceColor','r')
% hold on;
% plot(mean_params.grid,mean_params.prior_var,'linewidth',2,'Color','b')
% hold on;
% plot(mean_params.grid,mean_params.values,'linewidth',2,'Color','k')
% title([target_axis '-axis ' 'mean function'])
% xlabel('Stim location')
% 
% 
% % subplot(1,2,2)
% 
%  if output_plot
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 4 4];
% saveas(fig,[fig_name  'mean.png'])
% close(fig)
%  end
% 
% 
% %% Draw the predictive variance function:
% 
% fig=figure(5);
% 
% scatter(var_params.data.X,var_params.data.Y,...
%     'MarkerFaceColor','r')
% hold on;
% plot(var_params.grid,var_params.values,'linewidth',2,'Color','k')
% hold on;
% plot(var_params.grid, max(var_params.values)*(var_params.prior_var)/max(var_params.prior_var),'linewidth',2,'Color','b')
% 
% title([target_axis '-axis ' 'variance function'])
% xlabel('Stim location')
% % ylabel('Normalized current')
% 
% 
%  if output_plot
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 4 4];
% saveas(fig,[fig_name  'var.png'])
% close(fig)
%  end

end

