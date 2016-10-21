function psth = get_event_psths(events,t)


points = zeros(size(t));

for j = 1:length(events)
    these_events = events{j};
    if ~isempty(these_events)
        times = these_events(:,4);
        for i = 1:length(times)
            ind = find(t < times(i),1,'last')
            points(ind) = points(ind) + 1;
        end    
    end
end
points(end) = 0;
psth = smoothts(points,'b',100)*100;