function [] = plot_est_vs_truth(neurons,save_path)


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



figure_index=1; % can change to random numbers 



figure(figure_index)
scatter(gamma_truth,gamma_est,'Marker','o','SizeData',25,...
    'MarkerFaceColor','g', 'MarkerEdgeColor','g', 'MarkerFaceAlpha',0.8)
x=[0 1];y=[0 1];
hold on;
line(x,y,'Color','red','LineStyle','--')
hold off;

xlim([0 1]);
ylim([0 1]);

xlabel('True \gamma');
ylabel('Estimated \gamma');



figure(figure_index+1)
scatter(gain_truth(find(gamma_truth>0)),gain_est(find(gamma_truth>0)),'Marker','o','SizeData',25,...
    'MarkerFaceColor','g', 'MarkerEdgeColor','g', 'MarkerFaceAlpha',0.8)
x=[0 1];y=[0 1];
hold on;
line(x,y,'Color','red','LineStyle','--')
hold off;

xlim([0.00 0.03]);
ylim([0.00 0.03]);

xlabel('True optical gain');
ylabel('Estimated optical gain');



saveas(figure_index,strcat(save_path,'plots/', 'gamma_fits','.png'));
saveas(figure_index+1,strcat(save_path,'plots/', 'gain_fits','.png'));
close(figure_index,figure_index+1)



end