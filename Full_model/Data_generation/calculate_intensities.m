function [intensity] = calculate_intensities(stim,params,spike_curves,timepoints)
% Obtain the spike and psc event intensities
intensity=struct;
%intensity.spike = [];
%intensity.event=[];

[~, Ia]=min(abs(stim - spike_curves.current));
spike_param = struct;
spike_param.mean=spike_curves.mean(Ia);
spike_param.sd=spike_curves.sd(Ia);

% Spikes are generated from a truncated normal distribution, truncated at
% time_max
intensity.spike =  normpdf(timepoints,spike_param.mean,spike_param.sd);
intensity.event =  normpdf(timepoints,spike_param.mean+params.delay_mean, sqrt(spike_param.sd^2+params.delay_var) );
