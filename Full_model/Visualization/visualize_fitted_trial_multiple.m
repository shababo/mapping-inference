function [covered_flags] = visualize_fitted_trial_multiple(these_trials, yval,plot_params)
%%
% Visualize multiple trials at the same locations, ordered by power levels
if length(unique(yval))==1
   
plot_params.gap=yval(1)*0.2;
else
plot_params.gap=range(yval)/length(unique(yval));
end

n_trials = length(these_trials);
covered_flags=ones(length(unique( plot_params.loc_indices)),1);

if plot_params.by_neuron
   
    ytxt='Estimated stimulation size (log)';
else

    ytxt='Power level';
end
yval_max = max(yval);

% if plot_params.increasing 
% rankings_std= (max(yval)-yval)/range(yval);% scaled stim levels/rankings     
% else
% rankings_std= (yval-min(yval))/range(yval);% scaled stim levels/rankings 
% end
for i_trial=1:n_trials
    plot_params.shift= yval(i_trial);
    plot_params.this_loc = plot_params.these_indices(i_trial);
    this_trial = these_trials(i_trial);
    [covered_flag]=visualize_fitted_trial_single(this_trial, plot_params);
    if covered_flags(plot_params.this_loc)==1
        covered_flags(plot_params.this_loc)=covered_flag;
    end
end
if range(yval)>0
xlim([min(yval) max(yval)]+[-range(yval) range(yval)]*0.3);
else
xlim([min(yval) max(yval)]+[-max(yval) max(yval)]*0.3);    
end
ylabel('Time (ms)','FontSize', plot_params.lab_size); xlabel(ytxt,'FontSize', plot_params.lab_size);

