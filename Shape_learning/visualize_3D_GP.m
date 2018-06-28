function [ ]=visualize_3D_GP(GP_samples,plot_params)
%% Only draw the first sample! 
locations=GP_samples.full.locations;
threshold = plot_params.threshold;

unique_z=unique(locations(:,3));
n_grid = plot_params.n_grid;

%%
switch plot_params.type
    case 'scatter'
        this_sample=GP_samples.full.samples(:,1);
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
            % fig = gcf;
            % fig.PaperUnits = 'inches';
            % fig.PaperPosition = [0 0 4 4];
            % saveas(fig,[fig_name 'z'  num2str(unique_z(i_z)) '.png'])
            % close(fig)
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
            fig= figure(i_z)
            
            heatmap(Vq,x_grid,y_grid)
            caxis(cscale)
            xlabel('x');
            ylabel('y');
            title(['z plane:' num2str(unique_z(i_z)) ])
            
        end
end



