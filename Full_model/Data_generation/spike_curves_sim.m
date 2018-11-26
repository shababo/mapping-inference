function [spike_times,event_times] = spike_curves_sim(stim,params,spike_curves,varargin)
% ,spike_dist,event_times_dist <---- WE CAN RETURN THE WHOLE DISTRIBUTION
% HERE

if ~isempty(varargin) && ~isempty(varargin{1})
    time_max = varargin{1};
else
    time_max=140;%spike_curves.time_max;
end

spike_times = [];event_times=[];
[~, Ia]=min(abs(stim - spike_curves.current));
spike_param = struct;
spike_param.mean=spike_curves.mean(Ia);
spike_param.sd=spike_curves.sd(Ia);
spike_one = time_max +1;
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
