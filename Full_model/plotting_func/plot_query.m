function [] = plot_query(experiment_query,neighbourhoods,save_path)
%
% cell_shape=prior_info.template_cell.cell_shape; %can replace this with estimated shapes
%
% x_shape_grid=[-15 -10 -6 -3 0 3 6 10 15];
% y_shape_grid=[-15 -10 -6 -3 0 3 6 10 15];
% z_shape_grid=[-30 -20 -12 -6 0 6 12 20 30];

group_names = setdiff(fieldnames(experiment_query),'batch_ID');


% find the neighbourhood info from the query
for i_group = 1:length(group_names)
    if ~isempty(experiment_query.(group_names{i_group}))
        i_neighbourhood= experiment_query.(group_names{i_group}).neighbourhood_ID;
        break;
    end
end

% Summarize the query to create a table for trial counts & event counts for
% each unique locations 
query_summary=summarize_query(experiment_query);




figure_index=i_neighbourhood;
this_neighbourhood = neighbourhoods(i_neighbourhood);
cell_locations = reshape([this_neighbourhood.neurons(:).location],3,[])';
x=cell_locations(:,1);
y=cell_locations(:,2);
z=cell_locations(:,3);
xlims=[min(x) max(x)];
ylims=[min(y) max(y)];
zlims=[min(z) max(z)];


% for  trials:
markersize=30;
color_list=  [{'r'} {'k'}];

% for cells 
nucleisize=45;
secondary_color='g';
group_color_list=[{'m'} {'b'} {'y'} {'c'} ]; % these values should be defined in the get_structure()

%%
fig=figure(figure_index); 
for i_cell = 1:length(this_neighbourhood.neurons)
    this_cell = this_neighbourhood.neurons(i_cell);
        switch this_cell.group_ID
            case 'connected'
                markercolor=group_color_list{1};
                markeralpha=1;
            case 'undefined'
                markercolor=group_color_list{2};
                markeralpha=1;
            case 'disconnected'
                markercolor=group_color_list{3};
                markeralpha=1;
            case 'alive'
                markercolor=group_color_list{4};
                    markeralpha=1;
            case 'secondary'
                markercolor=secondary_color;
                markeralpha=0.2;
        end
    scatter3(this_cell.location(2),this_cell.location(1),this_cell.location(3),...
        'SizeData',nucleisize,...
        'MarkerEdgeColor',markercolor,...
        'MarkerFaceColor','w',...
        'MarkerFaceAlpha',markeralpha)
    
    hold on;
    
    
end


for i_group = 1:length(group_names)
    
    this_group=group_names{i_group};
   

    max_trial_count = max(query_summary.(this_group)(:,3));
    for i_loc = 1:size(query_summary.(this_group),1)
        
        loc_info=query_summary.(this_group)(i_loc,:);
        loc_coord = neighbourhoods(i_neighbourhood).neurons(loc_info(1)).stim_locations.(this_group).grid(loc_info(2),:);
        % for undefined group, bring the loc closer to the cell for easy
        % comparison 
        cell_coord = neighbourhoods(i_neighbourhood).neurons(loc_info(1)).location;
        if strcmp(this_group,'undefined')
            lift_flag=true;
           old_loc_coord = loc_coord;
           loc_coord(3)=cell_coord(3);
        else
           lift_flag=false;
        end
        
        if size(query_summary.(this_group),2)==4
            event_freq = loc_info(4)/loc_info(3);
            event_freq=min(1,event_freq);
        else
            event_freq = 1;
        end
        scatter3(loc_coord(2),loc_coord(1),loc_coord(3),...
                    'Marker','s',...
            'SizeData',markersize*(loc_info(3)/max_trial_count),...
            'MarkerEdgeColor',color_list{i_group},...
            'MarkerFaceColor',color_list{i_group},...
            'MarkerFaceAlpha',event_freq)
        hold on;
        if lift_flag
            line( [loc_coord(2) old_loc_coord(2)],[loc_coord(1) old_loc_coord(1)],[loc_coord(3) old_loc_coord(3)],...
            'LineStyle',':','LineWidth',1,...
            'Color',[0 0 0 0.2])
            hold on;
            scatter3(old_loc_coord(2),old_loc_coord(1),old_loc_coord(3),...
                'Marker','s',...
                'SizeData',markersize*(loc_info(3)/max_trial_count),...
                'MarkerEdgeColor',color_list{i_group},...
                'MarkerFaceColor','w',...
                'MarkerFaceAlpha',0,...
                'MarkerEdgeAlpha',0.2)
            
        end
    end
end
view(-20,20)
xlim(xlims);
ylim(ylims);
zlim(zlims);


for i_legend = 1:length(group_color_list)
  LegendHandels(i_legend) =  scatter3(nan,nan,nan,...
                'Marker','o',...
                'SizeData',nucleisize,...
                'MarkerEdgeColor',group_color_list{i_legend},...
                'MarkerFaceColor','w',...
                'MarkerEdgeAlpha',1);
            hold on;
end

LegendHandels(i_legend+1) =  scatter3(nan,nan,nan,...
    'Marker','s',...
    'SizeData',markersize,...
    'MarkerEdgeColor',color_list{1},...
    'MarkerFaceColor',color_list{1},...
    'MarkerEdgeAlpha',0.5);
hold on;


LegendHandels(i_legend+2) =  scatter3(nan,nan,nan,...
    'Marker','s',...
    'SizeData',markersize,...
    'MarkerEdgeColor',color_list{2},...
    'MarkerFaceColor',color_list{2},...
    'MarkerFaceAlpha',0.5);
hold on;


LegendHandels(i_legend+3) =  scatter3(nan,nan,nan,...
    'Marker','s',...
    'SizeData',markersize,...
    'MarkerEdgeColor',color_list{2},...
    'MarkerFaceColor','w',...
    'MarkerEdgeAlpha',0.2);
hold on;

LegendHandels(i_legend+4) =  scatter3(nan,nan,nan,...
    'Marker','s',...
    'SizeData',markersize,...
    'MarkerEdgeColor',color_list{2},...
    'MarkerFaceColor',color_list{2},...
    'MarkerFaceAlpha',0.2);
hold on;



LegendHandels(i_legend+5) =  scatter3(nan,nan,nan,...
    'Marker','s',...
    'SizeData',markersize,...
    'MarkerEdgeColor',color_list{2},...
    'MarkerFaceColor',color_list{2},...
    'MarkerFaceAlpha',0.8);
hold on;

legend(LegendHandels,{'connected cell', 'undefined cell','disconnected cell','alive cell',...
    'singlespot trials','multispot trials (lifted)','multispot trials (original)',...
    '20% responses','80% responses'});


fig.Units = 'inches';
fig.Position = [0 0 8 5];
hold off;
%%
saveas(figure_index,strcat(save_path,'plots/', 'TrialsBatch',num2str(experiment_query.batch_ID),'Neighbourhood',num2str(figure_index),'.png'));
close(figure_index)
end