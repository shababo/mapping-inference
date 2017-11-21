function [] = plot_one_neighbourhood(this_neighbourhood,save_path)

primary_color='b';
secondary_color='g';
primary_alpha=0.6;
secondary_alpha=0.3;

figure_index=this_neighbourhood.neighbourhood_ID; % can change to random numbers 

cell_locations = reshape([this_neighbourhood.neurons(:).location],3,[])';

primary_index=find([this_neighbourhood.neurons(:).primary_indicator]);
secondary_index = setdiff(1:length(this_neighbourhood.neurons), primary_index);



x=cell_locations(:,1);
y=cell_locations(:,2);
z=cell_locations(:,3);
xlims=[min(x) max(x)];
ylims=[min(y) max(y)];
zlims=[min(z) max(z)];

figure(figure_index)
ax1=axes('Position',[0 0 1 1],'Visible','off');
ax2=axes('Position',[0.1 0.2 .8 .8],'Visible','off');
axes(ax2)
scatter3(y(primary_index),x(primary_index),z(primary_index),...
    'MarkerEdgeColor',primary_color,...
    'MarkerFaceColor',primary_color,...
    'MarkerFaceAlpha',primary_alpha)

hold on;
scatter3(y(secondary_index),x(secondary_index),z(secondary_index),...
    'MarkerEdgeColor',secondary_color,...
    'MarkerFaceColor',secondary_color,...
    'MarkerFaceAlpha',secondary_alpha)

view(-30,10)
xlim(xlims);
ylim(ylims);
zlim(zlims);
axes(ax1)
title_text = {['Neighbourhood',' ',num2str(figure_index)]};
text(0.4,0.08,title_text,'fontsize',15)



saveas(figure_index,strcat(save_path,'plots/', 'neighbourhood',num2str(figure_index),'.png'));
close(figure_index)
end