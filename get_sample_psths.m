function psth = get_sample_psths(posteriors,bin_edges)


samples = [];

for j = 1:length(events)
    these_events = events{j};
    if ~isempty(these_events)
        times = these_events(:,4);
        for i = 1:length(times)
            ind = find(t < times(i),1,'last');
            points(ind) = 1;
        end    
    end
end
psth = smoothts(points,'g',300,100)*100;