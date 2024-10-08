% build structs to aggregate all spike + current data
clear results_current_all
clear results_spikes_tofit
clear results_current_all
count = 1;

power_cells = 1:length(result_current);
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

shape_cells = length(power_cells)+(1:length(result_xy));
for i_cell = 1:length(result_xy)
    
    results_current_all(count).these_powers = [result_xy(i_cell).these_x_power' result_xy(i_cell).these_y_power'];
    results_current_all(count).peak_current_means = [result_xy(i_cell).x_max_curr_means result_xy(i_cell).y_max_curr_means];
    results_spikes_all(count).these_powers = [result_xy(i_cell).these_x_power' result_xy(i_cell).these_y_power'];
    results_spikes_all(count).spike_time_means = [result_xy(i_cell).x_spike_time_means result_xy(i_cell).y_spike_time_means];
    results_spikes_all(count).spike_time_jitter = [result_xy(i_cell).x_spike_time_jitter result_xy(i_cell).y_spike_time_jitter];
    
    results_spikes_tofit(count).power = {results_spikes_all(count).these_powers'};
    results_spikes_tofit(count).spike_times = {result_xy(i_cell).spike_times'};
    results_spikes_tofit(count).spike_time_means = [result_xy(i_cell).x_spike_time_means result_xy(i_cell).y_spike_time_means];
    
    count = count + 1;
    
end

%%
% power_cells = 1:length(result_gain_depth);

for i_cell = 1:length(result_gain_depth)
%     results_current_all(count).these_powers = result_gain_depth(i_cell).spike_unique_powers;
%     results_current_all(count).peak_current_means = result_current(i_cell).peak_current_means;
%     results_spikes_all(count).these_powers = unique(result_spikes(i_cell).power{1});
%     results_spikes_all(count).spike_time_means = result_gain_depth(i_cell).spike_time_means;
%     results_spikes_all(count).spike_time_jitter = result_gain_depth(i_cell).spike_time_jitter;
    
    results_spikes_tofit(i_cell).power = {result_gain_depth(i_cell).spike_targ_power};
    results_spikes_tofit(i_cell).spike_times = {result_gain_depth(i_cell).spike_times};
    results_spikes_tofit(i_cell).spike_time_means = result_gain_depth(i_cell).spike_time_means;
    
%     count = count + 1;
end

%%
% power_cells = 1:28
% [spike_curves_power, x_current_power, y_spikes_power, spike_curve_cell_ids_power] = ...
%     get_spike_curves(results_current_all(power_cells),results_spikes_all(power_cells));
% 
% 
% [spike_curves_shape, x_current_shape, y_spikes_shape, spike_curve_cell_ids_shape] = ...
%     get_spike_curves(results_current_all(shape_cells),results_spikes_all(shape_cells));

[~, x_current_xy, y_spikes_xy]=get_spike_curves(results_current_all(29:end),results_spikes_all(1:29:end));

%%
% results_spikes needs fields: power (a cell of length 1 annoyingly),
% spike_times (same cell thing), [1

for this_cell = 1:length(result_gain_depth)
    this_cell
%     if any(this_cell == power_cells)
        [gain_mle_out(this_cell).gain, gain_mle_out(this_cell).gain_MC, gain_mle_out(this_cell).loglklh] = get_MLE_pilot(results_spikes_tofit(this_cell), spike_curves_power);
%     else
%         [gain_mle(this_cell), gain_MC, loglklh] = get_MLE_pilot(results_spikes_tofit(this_cell), spike_curves_power);
%     end
end




%%

% ni_mean=isnan(y_spike_mean) | y_spike_mean > 200 | x_current > 1200;
% xdata=x_current(~ni_mean);ydata=y_spike_mean(~ni_mean);
cells_to_plot = [1:35];
% colors = [repmat([0 0 1],28,1); repmat([0 1 0],7,1)];
% colors_scatter = colors(spike_curve_cell_ids,:);
h = figure;
subplot(121)
plot(spike_curves_power.current*1000,spike_curves_power.mean/20,'color',[0 0 1])
hold on

% plot(spike_curves_shape.current*1000,spike_curves_shape.mean/20,'color',[0 1 0])
hold on
scatter(x_current_power,y_spikes_power/20,[],[0 0 1])
hold on
scatter(x_current_shape,y_spikes_shape/20,[],[0 1 0])

% scatter(x_current_pow,y_spikes_pow/20)
% scatter(x_current_xy,y_spikes_xy/20)

ylim([0 10])
xlim([0 2500])
xlabel('mean peak current (pA)')
ylabel('mean spike time (msec)')
title('Fitting f(), 35 cells, 5 to 21 powers per cell')
subplot(122)
% figure
plot_spike_pred(results_spikes_tofit(cells_to_plot),spike_curves,gain_mle(cells_to_plot),colors)
ylim([0 15])
xlim([0 15])


%%

figure

scatter([result_spikes.fluor_val],gain_mle(1:28))
xlabel('Observed Fluor (a.u.)')
ylabel('Estimated Gain (a.u)')
title('Observed Fluor vs Estimated Gain')

