function [tr_stat] = calculate_timerescaling_stat(this_trial)
% Calculate the goodness-of-fit test statistics:
% Rules: 1) only count the first event 
%           2) use the max time if there is no event
tmp = cumsum(this_trial.fitted.intensity.event);
intensity_total=tmp(2,:);
if ~isempty(this_trial.event_times)
[~,im]=min(abs(this_trial.fitted.timepoints-this_trial.event_times(1)));
tr_stat=sum(intensity_total(1:im));
else
    tr_stat=sum(intensity_total);
end
