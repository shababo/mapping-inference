function [] = plot_neuron_locations(neurons,save_path)

figure_index=1; % can change to random numbers 
cell_locations = reshape([neurons(:).location],3,[])';


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
scatter3(y,x,z,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b',...
    'MarkerFaceAlpha',0.6)
view(-30,10)
xlim(xlims);
ylim(ylims);
zlim(zlims);
axes(ax1)
title_text = {'Neurons'};
text(0.4,0.08,title_text,'fontsize',15)



saveas(figure_index,strcat(save_path,'plots/', 'neurons','.png'));
close(figure_index)
end