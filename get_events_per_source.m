function dpmm_events_per_source = get_events_per_source(x) 

if isstruct(x)
%     dpmm_events_per_source = length(x(end).classes)/(x(end).num_classes-1); 
%     dpmm_events_per_source = max(x(end).counts);
    max_sum = 0;
    for i = 1:length(x)
        max_sum = max_sum + max(x(i).counts);
    end
    dpmm_events_per_source = max_sum/length(x);
else
    dpmm_events_per_source = 1; 
end
