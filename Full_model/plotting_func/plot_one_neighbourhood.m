function [] = plot_one_neighbourhood(this_neighbourhood,prior_info,save_path,varargin)

if ~isempty(varargin) 
    draw_shape = varargin{1};
else
    draw_shape = false;
end

if ~isempty(varargin) && length(varargin)>1
    draw_uncertainty = varargin{2};
else
    draw_uncertainty = false;
end
%%
group_names = unique({this_neighbourhood.neurons(:).group_ID});
% find the neighbourhood info from the query
i_neighbourhood=this_neighbourhood.neighbourhood_ID;
       
figure_index=i_neighbourhood;
cell_locations = reshape([this_neighbourhood.neurons(:).location],3,[])';
x=cell_locations(:,1);
y=cell_locations(:,2);
z=cell_locations(:,3);
xlims=[min(x) max(x)];
ylims=[min(y) max(y)];
zlims=[min(z) max(z)];

if draw_uncertainty
    % for uncertainty
    r_z=range(zlims)/20;
    r_x=range(xlims)/20;
end

if draw_shape 
    shape_grid = meshgrid(x_grid,y_grid,z_grid);
    stim_thresh=prior_info.induced_intensity.fire_stim_threshold;
    gain_bounds=struct;
    gain_bounds.low=group_profile.inference_params.bounds.gain(1);
    gain_bounds.up=group_profile.inference_params.bounds.gain(2);
    template_power=50;
    template_shape=prior_info.template_cell.cell_shape;
    
    x_grid = [-25 -15  -6 0 6  15 25];x_shift=41;
    y_grid=[-25 -15  -6 0 6  15 25];y_shift=41;
    z_grid=2*[-25 -15   -6  0 6  15 25];z_shift=91;
    pixelsize=8; % for drawing the cell shape (not implemented)

end

% for cells 
nucleisize=45;
secondary_color='g';
group_color_list=[{'m'} {'b'} {'y'} {'c'} ]; % these values should be defined in the get_structure()
maxpixelsize=8; % for drawing the cell shape (not implemented)

%% Obtain current estimates 

i_batch=this_neighbourhood.batch_ID;
neurons=this_neighbourhood.neurons;
properties={'PR_params'};summary_stat={'mean','lower_quantile','upper_quantile'};
temp_output=grab_values_from_neurons(i_batch,neurons,properties,summary_stat);

estimates=temp_output.PR_params;
%%
fig=figure(figure_index); 
for i_cell = 1:length(this_neighbourhood.neurons)
    this_cell = this_neighbourhood.neurons(i_cell);
    if ~this_cell.primary_indicator
        markercolor=secondary_color;
        markeralpha=0.1;
    else
        switch this_cell.group_ID
            case 'connected'
                markercolor=group_color_list{1};
            case 'undefined'
                markercolor=group_color_list{2};
            case 'disconnected'
                markercolor=group_color_list{3};
            case 'alive'
                markercolor=group_color_list{4};
        end
        markeralpha=estimates.mean(i_cell)*0.9+0.1;
    end
    scatter3(this_cell.location(2),this_cell.location(1),this_cell.location(3),...
        'SizeData',nucleisize,...
        'MarkerEdgeColor',markercolor,...
        'MarkerFaceColor',markercolor,...
        'MarkerFaceAlpha',markeralpha,...
        'MarkerEdgeAlpha',markeralpha)
    
    hold on;
    
    if draw_shape && strcmp(this_cell.group_ID,'connected')
        % we can change this to 
        alpha_gain=this_neighbourhood.neurons(i_cell).gain_params(end).alpha;
        beta_gain=this_neighbourhood.neurons(i_cell).gain_params(end).beta;
        temp=normrnd(alpha_gain,beta_gain,[group_profile.inference_params.MCsamples_for_posterior 1]);
        for i_x = x_grid
            for i_y = y_grid
                for i_z = z_grid
                     gain_thresh=stim_thresh/template_power/template_shape(i_x+x_shift,i_y+y_shift,i_z+z_shift);
                    if gain_thresh > gain_bounds.up
                       markeralpha=0; 
                    elseif gain_thresh < gain_bounds.low
                        markeralpha=1;
                    else
                        gain_trans= (gain_thresh-gain_bounds.low)/(gain_bounds.up-gain_bounds.low);
                        gain_logit = log(gain_trans/(1-gain_trans));
                        markeralpha = 1-normcdf(gain_logit,alpha_gain,beta_gain);
                    end
                    scatter3(this_cell.location(2)+i_y,this_cell.location(1)+i_x,this_cell.location(3)+i_z,...
                        'SizeData',pixelsize,...
                        'MarkerEdgeColor',markercolor,...
                        'MarkerFaceColor',markercolor,...
                        'MarkerFaceAlpha',markeralpha,...
                        'MarkerEdgeAlpha',markeralpha)
                    hold on;
                end
            end
        end
        
    end
 
    if draw_uncertainty
         	range_bar = [this_cell.location(3)+ 2*r_z*(estimates.lower_quantile(i_cell)-0.5), this_cell.location(3)+2*r_z*(estimates.upper_quantile(i_cell)-0.5)];
        line((this_cell.location(2))*ones(1,2),(this_cell.location(1)-r_x)*ones(1,2),range_bar,...
            'LineStyle','-','LineWidth',3,...
            'Color',[1 0 0 0.5])
        hold on;
        
        full_bar = [ [this_cell.location(3)- r_z this_cell.location(3)+ 2*r_z*(estimates.lower_quantile(i_cell)-0.5)]' ...
            [this_cell.location(3)+ r_z this_cell.location(3)+ 2*r_z*(estimates.upper_quantile(i_cell)-0.5)]'];
        line( (this_cell.location(2))*ones(2,2),(this_cell.location(1)-r_x)*ones(2,2),full_bar,...
            'LineStyle','-','LineWidth',3,...
            'Color',[0 0 0 0.2])
        hold on;
%         range(range_bar)
    end
    
    
end
view(-20,20)
xlim(xlims);
ylim(ylims);
zlim(zlims);

% 
% for i_legend = 1:length(group_color_list)
%   LegendHandels(i_legend) =  scatter3(nan,nan,nan,...
%                 'Marker','o',...
%                 'SizeData',nucleisize,...
%                 'MarkerEdgeColor',group_color_list{i_legend},...
%                 'MarkerFaceColor','w',...
%                 'MarkerEdgeAlpha',1);
%             hold on;
% end
% 
% LegendHandels(i_legend+1) =  scatter3(nan,nan,nan,...
%     'Marker','s',...
%     'SizeData',markersize,...
%     'MarkerEdgeColor',color_list{1},...
%     'MarkerFaceColor',color_list{1},...
%     'MarkerEdgeAlpha',0.5);
% hold on;
% 
% 
% LegendHandels(i_legend+2) =  scatter3(nan,nan,nan,...
%     'Marker','s',...
%     'SizeData',markersize,...
%     'MarkerEdgeColor',color_list{2},...
%     'MarkerFaceColor',color_list{2},...
%     'MarkerFaceAlpha',0.5);
% hold on;
% 
% 
% LegendHandels(i_legend+3) =  scatter3(nan,nan,nan,...
%     'Marker','s',...
%     'SizeData',markersize,...
%     'MarkerEdgeColor',color_list{2},...
%     'MarkerFaceColor','w',...
%     'MarkerEdgeAlpha',0.2);
% hold on;
% 
% LegendHandels(i_legend+4) =  scatter3(nan,nan,nan,...
%     'Marker','s',...
%     'SizeData',markersize,...
%     'MarkerEdgeColor',color_list{2},...
%     'MarkerFaceColor',color_list{2},...
%     'MarkerFaceAlpha',0.2);
% hold on;
% 
% 
% 
% LegendHandels(i_legend+5) =  scatter3(nan,nan,nan,...
%     'Marker','s',...
%     'SizeData',markersize,...
%     'MarkerEdgeColor',color_list{2},...
%     'MarkerFaceColor',color_list{2},...
%     'MarkerFaceAlpha',0.8);
% hold on;
% 
% legend(LegendHandels,{'connected cell', 'undefined cell','disconnected cell','alive cell',...
%     'singlespot trials','multispot trials (lifted)','multispot trials (original)',...
%     '20% responses','80% responses'});


fig.Units = 'inches';
fig.Position = [0 0 8 5];
hold off;
%%
saveas(figure_index,strcat(save_path,'plots/', 'neighbourhood',num2str(figure_index),'.png'));
close(figure_index)
end