function p_val = kstest_on_events(events,trace_length,offset)

num_events = sum(cellfun(@(x) size(x,1),events));
all_events = zeros(num_events,4);

if num_events > 1
    nulldist = makedist('Exponential','mu',trace_length/num_events); 

    events_combined = sort(all_events(:,4) - offset)';

    events_diff = events_combined(1);
    % size(events_combined)
    events_diff = [events_diff events_combined(2:end)-events_combined(1:end-1)];
    [h, p_val, k, c] = kstest(events_diff,'CDF',nulldist);
else
    p_val = 1.0;
end


