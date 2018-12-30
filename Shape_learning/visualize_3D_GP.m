function [ ]=visualize_3D_GP(GP_samples,plot_params)
%% Only draw the first sample!
locations=GP_samples.full.locations_unique;
threshold = plot_params.threshold;

unique_z=unique(locations(:,3));
n_grid = plot_params.n_grid;
fig_name = plot_params.fig_name;
%%
switch plot_params.type
    case 'scatter'
        this_sample=GP_samples.full.samples_unique(:,1);
        alpha_scale = ( (this_sample-min(this_sample)) ./range(this_sample)) ;
        alpha_scale(alpha_scale<threshold)=threshold;
        % fig_name=['./Figures/SimShape'];
        for i_z = 1:length(unique_z)
            this_z_indices = find(locations(:,3)==unique_z(i_z));
            xy_locations=locations(this_z_indices,1:2);
            fig= figure(i_z)
            for i = 1:length(this_z_indices)
                if alpha_scale(this_z_indices(i))>threshold
                    scatter(xy_locations(i,1),xy_locations(i,2),...
                        'MarkerFaceAlpha',alpha_scale(this_z_indices(i)),...
                        'MarkerFaceColor','k','MarkerEdgeColor','w')
                    hold on;
                end
                
            end
            xlim([min(xy_locations(:,1)) max(xy_locations(:,1))])
            ylim([min(xy_locations(:,2)) max(xy_locations(:,2))])
            title(['z plane:' num2str(unique_z(i_z)) ])
            xlabel('y')
            ylabel('x')
            
            %  if output_plot
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 4 4];
            saveas(fig,[fig_name 'z'  num2str(unique_z(i_z)) '.png'])
            close(fig)
            %  end
        end
        %% Alternatively: draw a heatmap with linear interpolation:
    case 'heatmap'
        cscale=   [min(GP_samples.full.samples(:,1)) max(GP_samples.full.samples(:,1))];
        
        
        for i_z = 1:length(unique_z)
            
            this_z_indices = find(locations(:,3)==unique_z(i_z));
            this_sample=GP_samples.full.samples(this_z_indices ,1);
            xy_locations=locations(this_z_indices,1:2);
            x_grid = (1:n_grid)*(range(xy_locations(:,1))/n_grid) +min(xy_locations(:,1));
            y_grid = (1:n_grid)*(range(xy_locations(:,2))/n_grid) +min(xy_locations(:,2));
            
            [Xq,Yq] = meshgrid(x_grid,y_grid);
            Vq = griddata(xy_locations(:,1),xy_locations(:,2),this_sample,Xq,Yq);
            fig= figure(i_z);
            
            imagesc(y_grid,x_grid,Vq',cscale)
            ylabel('x');
            xlabel('y');
            title(['z plane:' num2str(unique_z(i_z)) ])
            %  if output_plot
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 4 4];
            saveas(fig,[fig_name 'z'  num2str(unique_z(i_z)) '.png'])
            close(fig)
            
        end
    case 'contrast'
        
        cscale_shape=   [0 1];
        cscale_dev=   [-1 1];
        
        for i_z = 1:length(unique_z)
            
            this_z_indices = find(locations(:,3)==unique_z(i_z));
            this_sample=GP_samples.full.samples_unique(this_z_indices ,1);
            xy_locations=locations(this_z_indices,1:2);
            x_grid = (1:n_grid)*(range(xy_locations(:,1))/n_grid) +min(xy_locations(:,1));
            y_grid = (1:n_grid)*(range(xy_locations(:,2))/n_grid) +min(xy_locations(:,2));
            
            [Xq,Yq] = meshgrid(x_grid,y_grid);
            
            Vq = griddata(xy_locations(:,1),xy_locations(:,2),this_sample,Xq,Yq);
            Mq = griddata(xy_locations(:,1),xy_locations(:,2),GP_samples.full.mean(this_z_indices),Xq,Yq);
            
            
            %%
            fig= figure(i_z)
            
            subplot(1,3,1)
            imagesc(y_grid,x_grid,Vq',cscale_shape)
            ylabel('x (um)','FontSize',21);
            xlabel('y (um)','FontSize',21);
            ax = gca;
            ax.FontSize = 19;
            
            title(['Simulated cell'],'FontSize',21)
            %  if output_plot
            subplot(1,3,2)
            
            imagesc(y_grid,x_grid,Mq',cscale_shape)
%             ylabel('x (um)','FontSize',21);
            xlabel('y (um)','FontSize',21);
            ax = gca;
            ax.FontSize = 19;
             ax.YTick=[];
           
            title(['Mean'],'FontSize',21)
            hp3 = get(subplot(1,3,2),'Position');
            colorbar('Position', [hp3(1)+hp3(3)+0.01  hp3(2)  0.015  hp3(2)+hp3(3)*3.3])
            
            
            
            subplot(1,3,3)
            
            std_diag=sqrt(GP_samples.full.var);
            Sq = griddata(xy_locations(:,1),xy_locations(:,2),std_diag(this_z_indices),Xq,Yq);
            cscale_std = [min(std_diag) max(std_diag)];
            %               imagesc(x_grid,y_grid,(Vq-Mq)',cscale_shape)
            imagesc(y_grid,x_grid,Sq',cscale_std)
            %             ylabel('x');
%             ylabel('x (um)','FontSize',21);
            xlabel('y (um)','FontSize',21);
            ax = gca;
            ax.FontSize = 19;
            ax.YTick=[];
            title(['Standard Deviation'],'FontSize',21)
            hp4 = get(subplot(1,3,3),'Position');
            colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.015  hp4(2)+hp4(3)*3.3])
            
            % colorbar('Position', [hp4(1)+hp4(3)+0.1  hp4(2)  0.1  hp4(2)+hp4(3)*2.1])
            %%
            fig = gcf;
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 10 6];
            %             saveas(fig,[fig_name 'z'  num2str(unique_z(i_z)) '.png'])
            
            % print('-dpng','-r1500', [fig_name 'z'  num2str(unique_z(i_z))])
            %             close(fig)
            
        end
end



