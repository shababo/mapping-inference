function [] = visualize_fitted_trial_multiple(these_trials, rankings,plot_params)

% Visualize multiple trials at the same locations, ordered by power levels
vertical_gap=plot_params.vertical_gap;
n_trials = length(these_trials);
if plot_params.by_neuron
   
    ytxt='Trials ordered by est. stim.';
else

    ytxt='Trials ordered by power';
end
rankings_std= (rankings-min(rankings))/range(rankings);% scaled stim levels/rankings 
n_ranks = length(unique(rankings_std));
for i_trial=1:n_trials
    plot_params.vertical_shift= vertical_gap*n_ranks*(1-rankings_std(i_trial));
    this_trial = these_trials(i_trial);
    visualize_fitted_trial_single(this_trial, plot_params)
end
ylim([0  (n_ranks+1)*vertical_gap]);
xlabel('Time (ms)'); ylabel(ytxt);

