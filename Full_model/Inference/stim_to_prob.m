function [first_spike_probability] = stim_to_prob(...
    stim,stim_unique,prob_trace)
    [~, stim_idx] = min(abs(stim-stim_unique));
    first_spike_probability = prob_trace(stim_idx);
end