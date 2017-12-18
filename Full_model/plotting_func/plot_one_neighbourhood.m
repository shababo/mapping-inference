function [] = plot_one_neighbourhood(neighbourhood,varargin)

if ~isempty(varargin) && ~isempty(varargin{1})
    fighandle = varargin{1};
else
    fighandle = gcf;
end

%% Draw three plots
% Max projection, 3D location with group colors, Uncertainty of gain and
% gamma 
default_color={'k'};
default_fill_color={'w'};

group_colors={'r' 'g' 'c' 'b' 'k'}; % disconnected, undefined, connected, alive, secondary
fill_colors={'r' 'g' 'c' 'b' 'k'}; % disconnected, undefined, connected, alive, secondary
group_names={'disconnected','undefined','connected','alive','secondary'};
neuron_alpha=0.5; % transparency
neuron_size=20;

    barcolor=[0 0 0 0.2];
    r=15; %length of bars
    barwidth=3;
    bargap=2;
    gain_bounds=[0.005 0.03];
    
fig_size=[0 0 15 5];
%


neurons=neighbourhood.neurons;
cell_locations = reshape([neurons(:).location],3,[])';
x=cell_locations(:,1);
y=cell_locations(:,2);
z=cell_locations(:,3);
xlims=[-152 152];
ylims=[-152 152];
zlims=[min(cell_locations(:,3)) max(cell_locations(:,3))];

i_batch=neighbourhood.batch_ID;
properties={'PR_params','gain_params'};summary_stat={'lower_quantile','mean','upper_quantile'};
posteriors=grab_values_from_neurons(i_batch,neurons,properties,summary_stat);


%% Draw 2D max proj

figure(fighandle);

subplot(3,4,(neighbourhood.neighbourhood_ID - 1)*4 + 1) 
if isfield(neighbourhood,'stack')
    imagesc(max(neighbourhood.stack,[],3));
end
title(['Neighbourhood ' num2str(neighbourhood.neighbourhood_ID)])
%%
% subplot(3,4,(neighbourhood.neighbourhood_ID - 1)*4 + 2) 
% 
% for i_cell = 1:length(neurons)
%     switch neighbourhood.neurons(i_cell).group_ID{end}
%         case 'disconnected'
%            i_group = 1;
%         case 'undefined'
%             i_group = 2;
%         case 'connected'
%             i_group = 3;
%         case 'alive'
%             i_group = 4;
%         case 'secondary'
%             i_group = 5;
%     end
%      edgecolor=group_colors{i_group};
%             fillcolor=fill_colors{i_group};
% scatter(x(i_cell),y(i_cell),neuron_size,edgecolor,'filled');
% %     'MarkerFaceColor',fillcolor);%,...
% %     'MarkerFaceAlpha',neuron_alpha)
% hold on;
% end
% xlim(xlims);
% ylim(ylims);
% xlabel('x (um)')
% ylabel('y (um)')
% 
% hold off;

%%
% Draw 3D loc 
subplot(3,4,(neighbourhood.neighbourhood_ID - 1)*4 + 3) 

for i_cell = 1:length(neurons)
    switch neighbourhood.neurons(i_cell).group_ID{end}
        case 'disconnected'
           i_group = 1;
        case 'undefined'
            i_group = 2;
        case 'connected'
            i_group = 3;
        case 'alive'
            i_group = 4;
        case 'secondary'
            i_group = 5;
    end
     edgecolor=group_colors{i_group};
            fillcolor=edgecolor;
scatter3(x(i_cell),y(i_cell),z(i_cell),neuron_size,edgecolor,'filled');
% scatter3(x(i_cell),y(i_cell),z(i_cell),...
%     'SizeData',neuron_size,...
%     'MarkerEdgeColor',edgecolor,...
%     'MarkerFaceColor',fillcolor);%,...
%     'MarkerFaceAlpha',neuron_alpha)
hold on;
end
view(-30,10)
xlim(xlims);
ylim(ylims);
zlim(zlims);
xlabel('x (um)')
ylabel('y (um)')
zlabel('z (um)')

hold off;

% Draw uncertainty 
%%
% subplot(3,4,(neighbourhood.neighbourhood_ID - 1)*4 + 4)
% 
% for i_cell = 1:length(neurons)
%     switch neighbourhood.neurons(i_cell).group_ID{end}
%         case 'disconnected'
%            i_group = 1;
%         case 'undefined'
%             i_group = 2;
%         case 'connected'
%             i_group = 3;
%         case 'alive'
%             i_group = 4;
%         case 'secondary'
%             i_group = 5;
%     end
%     
%      edgecolor=group_colors{i_group};
%      fillcolor=default_fill_color{1};
%     
% 
% xloc=x(i_cell)*ones(2,1)-barwidth/2-bargap;
% yloc =[r*(-0.5),  r*(0.5)]+y(i_cell);
% line(xloc,yloc,...
%     'LineStyle','-','LineWidth',barwidth,...
%     'Color',barcolor)
% hold on;
% xloc=x(i_cell)*ones(2,1)-barwidth/2-bargap;
% yloc =[r*( posteriors.PR_params.lower_quantile(i_cell)-0.5),  r*(  posteriors.PR_params.upper_quantile(i_cell)-0.5)]+...
%     y(i_cell);
% line(xloc,yloc,...
%     'LineStyle','-','LineWidth',barwidth,...
%     'Color',edgecolor)
% hold on;
% 
% xloc=x(i_cell)*ones(2,1)+barwidth/2+bargap;
% yloc =[r*(-0.5),  r*(0.5)]+y(i_cell);
% line(xloc,yloc,...
%     'LineStyle','-','LineWidth',barwidth,...
%     'Color',barcolor)
% hold on;
% xloc=x(i_cell)*ones(2,1)+barwidth/2+bargap;
% yloc =[r* (( posteriors.gain_params.lower_quantile(i_cell)-gain_bounds(1))/range(gain_bounds)-0.5),  ...
%     r*((gain_bounds(2)- posteriors.gain_params.lower_quantile(i_cell))/range(gain_bounds)-0.5)]+...
%     y(i_cell);
% line(xloc,yloc,...
%     'LineStyle','-','LineWidth',barwidth,...
%     'Color',edgecolor)
% hold on;
% 
% end
% 
% 
% 
% % for i_legend = 1:length(group_colors)
% %   LegendHandels(i_legend) =  scatter(nan,nan,...
% %                 'Marker','o',...
% %                 'SizeData',neuron_size,...
% %                 'MarkerEdgeColor',group_colors{i_legend},...
% %                 'MarkerFaceColor',group_colors{i_legend});%,...
% % %                 'MarkerEdgeAlpha',neuron_alpha);
% %             hold on;
% % end
% % legend(LegendHandels,group_names,'Location','southeast');
% 
% xlim(xlims);
% ylim(ylims);
% xlabel('x (um)')
% ylabel('y (um)')
% % hold off;

%%
fig.Units = 'inches';
fig.Position =fig_size;

%%
% saveas(figure_index,strcat('plots/', 'Summary_n',num2str(figure_index),'_b',num2str(i_batch),'.png'));
% close(figure_index)
