function [] = plot_est_vs_truth(this_neighbourhood,save_path,varargin)

neurons = this_neighbourhood(:).neurons;


if ~isempty(varargin)
   figure_index = varargin{1}; 
else
    figure_index=1; % can change to random numbers 
end
number_of_cells = length(neurons);
gain_truth=zeros(number_of_cells,1);
gamma_truth=zeros(number_of_cells,1);
gain_est=zeros(number_of_cells,1);
gamma_est=zeros(number_of_cells,1);

for i_cell = 1:number_of_cells
   gain_truth(i_cell)=neurons(i_cell).truth.optical_gain; 
   gamma_truth(i_cell)=neurons(i_cell).truth.PR;
   gain_est(i_cell)=neurons(i_cell).gain_params(end).mean;
   gamma_est(i_cell)=neurons(i_cell).PR_params(end).mean;
end



group_color_list=[{'g'} {'b'} {'r'} {'c'}]; % these values 

fig=figure(figure_index);
subplot(1,2,1)       % add first plot in 2 x 1 grid
for i=1:length(neurons)
    switch neurons(i).group_ID
        case 'undefined'
            
            markercolor=group_color_list{1};
        case 'connected'
            
            markercolor=group_color_list{2};
        case 'disconnected'
            
            markercolor=group_color_list{3};
        case 'alive'
            markercolor=group_color_list{4};
    end
scatter(gamma_truth(i),gamma_est(i),'Marker','o','SizeData',25,...
    'MarkerFaceColor',markercolor, 'MarkerEdgeColor',markercolor, 'MarkerFaceAlpha',0.8)
hold on;
end
x=[0 1];y=[0 1];

line(x,y,'Color','red','LineStyle','--')
hold off;

xlim([0 1]);
ylim([0 1]);

xlabel('True \gamma');
ylabel('Estimated \gamma');


subplot(1,2,2)       % add first plot in 2 x 1 grid


scatter(gain_truth(find(gamma_truth>0)),gain_est(find(gamma_truth>0)),'Marker','o','SizeData',25,...
    'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.8)
hold on;
x=[0 1];y=[0 1];
line(x,y,'Color','red','LineStyle','--')
hold off;

xlim([0.00 0.03]);
ylim([0.00 0.03]);

xlabel('True optical gain');
ylabel('Estimated optical gain');



fig.Units = 'inches';
fig.Position = [0 0 8 4];

saveas(figure_index,strcat(save_path,'plots/', 'Fits',num2str(figure_index),'.png'));
close(figure_index)



end