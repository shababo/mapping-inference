function [pred_spike_time]=predict_spike_time(stim,spike_curves)

[~, iloc]=min( abs(stim-spike_curves.current));
pred_spike_time=spike_curves.mean(iloc);