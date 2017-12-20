function [figure_handle] = plot_intensity_one_trial(trial_intensity,figure_handle, varargin)

% varargin for color and other graph specs
color_map = cool(length(trial_intensity.stim_neurons)+1);
xgrid = (1:length(trial_intensity.estimated_intensity))/20;

%
figure(figure_handle)
plot(xgrid,trial_intensity.estimated_intensity,'LineWidth',3,'Color','k')
hold on;
ymax=max(trial_intensity.estimated_intensity);
if ~isempty(trial_intensity.event_times)
    for i_event = 1:length(trial_intensity.event_times)
        line(xgrid([trial_intensity.event_times(i_event),trial_intensity.event_times(i_event)]),[0,ymax],...
            'LineStyle',':','LineWidth',3,'Color','r')
        hold on;
    end
end
for i_cell = 1:length(trial_intensity.stim_neurons)
    temp_intensity=trial_intensity.stim_neurons(i_cell).intensity*trial_intensity.stim_neurons(i_cell).PR;
    plot(xgrid,temp_intensity,'LineWidth',3,'Color',color_map(i_cell,:))
    hold on;
end
 plot(xgrid,trial_intensity.background_rate*ones(length(trial_intensity.estimated_intensity),1),'LineWidth',3,'Color',color_map(end,:))
    hold on;

    LegendHandels(1)= plot(nan,nan,'LineWidth',3,'Color','k');
    hold on;
    legend_names{1}='Overall';
    LegendHandels(2)= plot(nan,nan,'LineStyle',':','LineWidth',3,'Color','r');
    hold on;
    legend_names{2}='Events';
    

    for i_legend = 1:length(trial_intensity.stim_neurons)
        LegendHandels(i_legend+2) =  ...
            plot(nan,nan,'LineWidth',3,'Color',color_map(i_legend,:));
        hold on;
        legend_names{i_legend+2}=['Neuron', ' ', num2str(trial_intensity.stim_neurons(i_legend).cell_ID)];
    end
    LegendHandels(i_legend+3) =  ...
        plot(nan,nan,'LineWidth',3,'Color',color_map(end,:));
    hold on;
    legend_names{i_legend+3}='Background';
    legend(LegendHandels,legend_names,'Location','northeast');

    xlabel('Time (ms)', 'FontSize',14)
    ylabel('Intensity', 'FontSize',14)
    title(['Trial', ' ',num2str(trial_intensity.trial_ID)],'FontSize',14)
    
