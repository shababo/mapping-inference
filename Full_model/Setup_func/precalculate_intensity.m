
function induced_intensity= precalculate_intensity(induced_intensity,...
    template_cell,delay,current_template)


[induced_intensity.intensity_grid,~] = get_first_spike_intensity(...
    current_template,induced_intensity.stim_grid,...
    template_cell,delay);
induced_intensity.probability_grid=sum(induced_intensity.intensity_grid,2);

induced_intensity.minimum_stim_threshold=induced_intensity.stim_grid(min(find(induced_intensity.probability_grid>0.01)));
induced_intensity.fire_stim_threshold=induced_intensity.stim_grid(min(find(induced_intensity.probability_grid>0.99)));



% %%
% induced_intensity=experiment_setup.prior_info.induced_intensity;
% template_cell=  experiment_setup.prior_info.template_cell;
% delay=experiment_setup.prior_info.delay;
% current_template=experiment_setup.prior_info.current_template;