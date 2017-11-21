function [] = plot_designs(experiment_query,neighbourhoods,save_path)
%
% cell_shape=prior_info.template_cell.cell_shape; %can replace this with estimated shapes
%
% x_shape_grid=[-15 -10 -6 -3 0 3 6 10 15];
% y_shape_grid=[-15 -10 -6 -3 0 3 6 10 15];
% z_shape_grid=[-30 -20 -12 -6 0 6 12 20 30];
%
%

group_names = setdiff(fieldnames(experiment_query),'batch_ID');


primary_color='b'; secondary_color='g';
primary_alpha=0.6; secondary_alpha=0.3;

% find the neighbourhood info from the query
for i_group = 1:length(group_names)
    if ~isempty(experiment_query.(group_names{i_group}))
        i_neighbourhood= experiment_query.(group_names{i_group}).neighbourhood_ID;
        break;
    end
end

figure_index=i_neighbourhood;
this_neighbourhood = neighbourhoods(i_neighbourhood);
cell_locations = reshape([this_neighbourhood.neurons(:).location],3,[])';



primary_index=find([this_neighbourhood.neurons(:).primary_indicator]);
secondary_index = setdiff(1:length(this_neighbourhood.neurons), primary_index);




x=cell_locations(:,1);
y=cell_locations(:,2);
z=cell_locations(:,3);
xlims=[min(x) max(x)];
ylims=[min(y) max(y)];
zlims=[min(z) max(z)];

markershape='o';
markeralpha=0.5;
markersize=10;
color_list=  [{'r'} {'k'} {'m'} {'y'} {'c'}];

nucleisize=20;
maxpixelsize=8;

%
figure(figure_index)
ax1=axes('Position',[0 0 1 1],'Visible','off');
ax2=axes('Position',[0.1 0.2 .8 .8],'Visible','off');
axes(ax2)

for i_group = 1:length(group_names)
    
    cells_this_group= find(get_group_inds(this_neighbourhood,group_names{i_group}));
    query_this_group=experiment_query.(group_names{i_group});
    
    
    primary_this_group = intersect(primary_index,cells_this_group);
    scatter3(y(primary_this_group),x(primary_this_group),z(primary_this_group),...
        'SizeData',nucleisize,...
        'MarkerEdgeColor',primary_color,...
        'MarkerFaceColor',primary_color,...
        'MarkerFaceAlpha',primary_alpha)
    
    hold on;
    
    secondary_this_group = intersect(secondary_index,cells_this_group);
    scatter3(y(secondary_this_group),x(secondary_this_group),z(secondary_this_group),...
        'SizeData',nucleisize,...
        'MarkerEdgeColor',secondary_color,...
        'MarkerFaceColor',secondary_color,...
        'MarkerFaceAlpha',secondary_alpha)
    
    hold on;
    
    
    for i_trial = 1:length(query_this_group.trials)
        x_trial =query_this_group.trials(i_trial).locations(:,1);
        y_trial =query_this_group.trials(i_trial).locations(:,2);
        z_trial =query_this_group.trials(i_trial).locations(:,3);
        
        scatter3(y_trial,x_trial,z_trial,...
            'SizeData',markersize,...
            'MarkerEdgeColor',color_list{i_group},...
            'MarkerFaceColor',color_list{i_group},...
            'MarkerFaceAlpha',markeralpha)
        hold on;
    end
    
end

view(-30,10)
xlim(xlims);
ylim(ylims);
zlim(zlims);
axes(ax1)
title_text = {'Stimulation locations'};
axes(ax1)
text(0.2,0.08,title_text,'fontsize',15)


%%
saveas(figure_index,strcat(save_path,'plots/', 'TrialsBatch',num2str(experiment_query.batch_ID),'Neighbourhood',num2str(figure_index),'.png'));
close(figure_index)
end