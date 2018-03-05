% build structs to aggregate all spike + current data
clear results_current_all
clear results_spikes_tofit
clear results_current_all
count = 1;

for i_cell = 1:length(result_current)
    results_current_all(count).these_powers = unique(result_current(i_cell).power{1});
    results_current_all(count).peak_current_means = result_current(i_cell).peak_current_means;
    results_spikes_all(count).these_powers = unique(result_spikes(i_cell).power{1});
    results_spikes_all(count).spike_time_means = result_spikes(i_cell).spike_time_means;
    results_spikes_all(count).spike_time_jitter = result_spikes(i_cell).spike_time_jitter;
    
    results_spikes_tofit(count).power = result_spikes(i_cell).power;
    results_spikes_tofit(count).spike_times = result_spikes(i_cell).spike_times;
    results_spikes_tofit(count).spike_time_means = result_spikes(i_cell).spike_time_means;
    
    count = count + 1;
end

for i_cell = 1:length(result_xy)
    
    results_current_all(count).these_powers = [result_xy(i_cell).these_x_power' result_xy(i_cell).these_y_power'];
    results_current_all(count).peak_current_means = [result_xy(i_cell).x_max_curr_means result_xy(i_cell).y_max_curr_means];
    results_spikes_all(count).these_powers = [result_xy(i_cell).these_x_power' result_xy(i_cell).these_y_power'];
    results_spikes_all(count).spike_time_means = [result_xy(i_cell).x_spike_time_means result_xy(i_cell).y_spike_time_means];
    results_spikes_all(count).spike_time_jitter = [result_xy(i_cell).x_spike_time_jitter result_xy(i_cell).y_spike_time_jitter];
    
    results_spikes_tofit(count).power = {result_xy(i_cell).spatial_adj_power};
    results_spikes_tofit(count).spike_times = {result_xy(i_cell).spike_times};
    results_spikes_tofit(count).spike_time_means = [result_xy(i_cell).x_spike_time_means result_xy(i_cell).y_spike_time_means];
    
    count = count + 1;
    
end
%%

[spike_curves, x_current, y_spikes]=get_spike_curves(results_current_all,results_spikes_all);
%%
% results_spikes needs fields: power (a cell of length 1 annoyingly),
% spike_times (same cell thing), [1

for this_cell = 1:length(results_spikes_tofit)
    this_cell
    [gain_mle(this_cell), gain_MC, loglklh]=get_MLE_pilot(results_spikes_tofit(this_cell), spike_curves);
end




%%

% ni_mean=isnan(y_spike_mean) | y_spike_mean > 200 | x_current > 1200;
% xdata=x_current(~ni_mean);ydata=y_spike_mean(~ni_mean);
cells_to_plot = [1:28];
h = figure;
subplot(121)
plot(spike_curves.current*1000,spike_curves.mean/20)
hold on
scatter(x_current,y_spikes/20)
ylim([0 10])
xlim([0 2500])
xlabel('mean peak current (pA)')
ylabel('mean spike time (msec)')
title('Fitting f(), 35 cells, 5 to 21 powers per cell')
subplot(122)
% figure
colors = lines(35);%[repmat([0 0 1],28,1); repmat([0 1 0],7,1)];
plot_spike_pred(results_spikes_tofit(cells_to_plot),spike_curves,gain_mle(cells_to_plot),colors(cells_to_plot,:))
ylim([0 15])
xlim([0 15])


%%

figure

scatter([result_spikes.fluor_val],gain_mle(1:28))
xlabel('Observed Fluor (a.u.)')
ylabel('Estimated Gain (a.u)')
title('Observed Fluor vs Estimated Gain')

