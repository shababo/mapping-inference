
function induced_intensity= precalculate_intensity(induced_intensity,...
    template_cell,delay,current_template)


[induced_intensity.intensity_grid,~] = get_first_spike_intensity(...
    current_template,induced_intensity.stim_grid,...
    template_cell,delay);
induced_intensity.probility_grid=sum(induced_intensity.intensity_grid,2);

induced_intensity.minimum_stim_threshold=induced_intensity.stim_grid(min(find(induced_intensity.probility_grid>0.01)));
induced_intensity.fire_stim_threshold=induced_intensity.stim_grid(min(find(induced_intensity.probility_grid>0.99)));