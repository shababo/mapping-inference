function [] = visualize_fitted_trial_single(this_trial, plot_params)
alpha=plot_params.alpha;
if plot_params.prediction
    chosen_field = 'predicted';
else
      chosen_field = 'fitted';
end


realtimepoints=this_trial.(chosen_field).timepoints/this_trial.(chosen_field).time_factor;
n_neuron= size(this_trial.(chosen_field).intensity.event,1); % This will always include the background rate 
neuron_colors = lines(n_neuron);


if plot_params.spike 
    fits = this_trial.(chosen_field).intensity.spike;
else
    fits =this_trial.(chosen_field).intensity.event;
end

if plot_params.stack 
    % Stack the intensities 
    fits = cumsum(fits);
end
if isfield(plot_params, 'LineWidth')
    lw=plot_params.LineWidth;
else
    lw=5;
end
if isfield(plot_params,'vertical_shift')
    vertical_shift=plot_params.vertical_shift;
else
    vertical_shift=0;
end 
if    isfield(plot_params,'vertical_gap') 
    barheight=plot_params.vertical_gap;
else
    barheight = max(max(fits));
end

for i = 1:n_neuron
    plot(realtimepoints, fits(i,:)+vertical_shift*ones(1,size(fits,2)),'Color',[neuron_colors(i,:) alpha],'LineWidth',lw)
    hold on;
end
if isfield(this_trial,'event_times')
if ~isempty(this_trial.event_times)
tmp=this_trial.event_times/this_trial.(chosen_field).time_factor;
plot([tmp;tmp], [vertical_shift;vertical_shift+barheight],'-','LineWidth',lw+1,'Color','k')
end
end

if ~isfield(plot_params,'vertical_shift')
    xlabel('Time (ms)');
    if plot_params.prediction
    ylabel('Predicted intensity')
    else
        ylabel('Fitted intensity')
    end
end

   
% Visualize multiple trials at the same locations, ordered by power levels 