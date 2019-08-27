function [] = visualize_fitted_trial_multiple(these_trials, rankings,plot_params)

% Visualize multiple trials at the same locations, ordered by power levels
vertical_gap=plot_params.vertical_gap;
n_trials = length(these_trials);
if plot_params.by_neuron
   
    ytxt='Trials ordered by est. stim.';
else

    ytxt='Trials ordered by power';
end
for i_trial=1:n_trials
    plot_params.vertical_shift= vertical_gap*(n_trials - rankings(i_trial));
    this_trial = these_trials(i_trial);
    visualize_fitted_trial_single(this_trial, plot_params)
end
ylim([0  (n_trials)*vertical_gap]);
xlabel('Time (ms)'); ylabel(ytxt);

