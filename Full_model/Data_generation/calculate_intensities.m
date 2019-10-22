function [intensity] = calculate_intensities(stim,delay_params,spike_curves,timepoints,params)
% Obtain the spike and psc event intensities
intensity=struct;
%intensity.spike = [];
%intensity.event=[];
spike_param = struct;
spike_param.mean=spike_curves.F_mean(spike_curves.mean_param,stim/spike_curves.current_multiplier);
    spike_param.sd= spike_curves.inflate_func(stim, params.inflate_func_coef);

% Spikes are generated from a truncated normal distribution, truncated at
% time_max
intensity.spike =  normpdf(timepoints,spike_param.mean,spike_param.sd);
intensity.event =  normpdf(timepoints,spike_param.mean+delay_params.delay_mean, sqrt(spike_param.sd^2+delay_params.delay_var) );
