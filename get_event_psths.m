function psth = get_event_psths(events,t,amp_bounds)


points = zeros(size(t));

for j = 1:length(events)
    these_events = events{j};
    if ~isempty(these_events)
        times = these_events(:,4);
        for i = 1:length(times)
            if times(i) > t(1) && times(i) < times(end) ...
                    && these_events(i,1) > amp_bounds(1) && ...
                    these_events(i,1) < amp_bounds(2)
                ind = find(t < times(i),1,'last')
                points(ind) = points(ind) + 1;
            end
        end    
    end
end
points(end) = 0;
psth = smoothts(points,'g',200,50)*100;