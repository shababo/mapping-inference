function [] = visualize_fitted_heatmap(these_trials, plot_params)
markerSize=plot_params.markerSize;
event_range=plot_params.event_range;
n_trials =length(these_trials);

first_event_time = zeros(n_trials,1);

for i_trial = 1:n_trials
tmp=cumsum(these_trials(i_trial).predicted.intensity.event);
intensity_total = tmp(end,:);
[~,im]=max(intensity_total);
first_event_time(i_trial)=these_trials(i_trial).predicted.timepoints(im)/these_trials(i_trial).predicted.time_factor;
end

[power_unique, ~, ind_power] =unique([these_trials(:).power_levels]);
for i =  1:length(power_unique)
    tmp_trials=these_trials(ind_power==i);
    tmp_time = first_event_time(ind_power==i);
    tmp_loc=reshape([tmp_trials.locations],3,[])';
    [locations_unique, ia, ic] = unique(tmp_loc,'rows');
    tmp_time_unique=tmp_time(ia);
    figure(i)
    scatter3(locations_unique(:,1),locations_unique(:,2),locations_unique(:,3),markerSize,tmp_time_unique,'filled')
    oldcmap = colormap;
    colormap( flipud(oldcmap) );
    colorbar;
    caxis(plot_params.event_range)
end
