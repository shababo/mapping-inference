function [spike_times,event_times] = spike_curves_sim(stim,params,spike_curves,varargin)
% Draw new spike and new events 

if ~isempty(varargin) && ~isempty(varargin{1})
    time_max = varargin{1};
else
    time_max=spike_curves.time_max;
end

spike_times = [];event_times=[];

stim=stim/spike_curves.current_multiplier;
spike_param = struct;
spike_param.mean=spike_curves.F_mean(spike_curves.mean_param,stim);
if params.combined_variance
    spike_param.sd= spike_curves.inflate_func(...
        sqrt(spike_curves.F_sd(spike_curves.sd_param,stim)^2+spike_curves.F_dev(spike_curves.dev_param,stim)^2), params.inflate_func_coef);
else
    spike_param.sd=spike_curves.inflate_func(...
        spike_curves.F_sd(spike_curves.sd_param,stim), params.inflate_func_coef);
end
% Spikes are generated from a truncated normal distribution, truncated at
% time_max
spike_one = normrnd(spike_param.mean,spike_param.sd);
% spike_one = normrnd(spike_param.mean,spike_param.sd);
delay_one = normrnd(params.delay_mean,sqrt(params.delay_var));
if spike_one < time_max
    spike_times = [spike_times spike_one];
    event_times = [event_times spike_one+delay_one];
else
    spike_times = [spike_times Inf];
    event_times = [event_times Inf];    
end
