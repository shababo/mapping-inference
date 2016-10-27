function dpmm_out = get_source_clusters(events)


num_events = sum(cellfun(@(x) size(x,1),events));
all_events = zeros(num_events,4);

event_i = 1;
for j = 1:length(events)

    if ~isempty(events{j})
        all_events(event_i:event_i+size(events{j},1)-1,:) = events{j};
    end
    event_i = event_i+size(events{j},1);
    
    
end

size(all_events)
assignin('base','loc_events',all_events)

all_events = all_events(:,[4]);
if num_events > size(all_events,2)
%     [dpmm_params, dpmm_gammas, dpmm_assign] = dpmm(all_events,25);
    dpmm_out = dpmm(all_events,25);
else
    dpmm_out = num_events;
    dpmm_gammas = [];
    dpmm_assign = [];
end