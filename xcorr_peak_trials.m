function [xcorr_peaks, test_stat] = xcorr_peak_trials(traces1,traces2,t_bounds,null)

if ~isempty(t_bounds)
    if null
        order = randperm(size(traces1,1));
    else
        order = 1:size(traces1,1);
    end
    traces1 = traces1(order,t_bounds(1):t_bounds(2));
    traces2 = traces2(:,t_bounds(1):t_bounds(2));
end




[this_xcorr,lags] = xcorr(traces1(:),traces2(:));
% figure; plot(lags,this_xcorr)
this_xcorr = this_xcorr(lags > -20 & lags < 20);
[xcorr_peaks, i] = max(this_xcorr);
assignin('base','lags',lags)
num_lags = length(this_xcorr);
test_stat = this_xcorr.*sqrt((num_lags - 2)./(1 - this_xcorr.^2));
test_stat = test_stat(i);

    

