params.g = .1;
% params.current_template = template;
clear responses
clear stims
loc_ind = 1;
spike_times = cells(1).spike_data(loc_ind).spike_times;
powers = cells(1).spike_data(loc_ind).powers;
count = 1;
for i = 1:length(spike_times)
    
    for j = 1:length(spike_times{i})
        responses(count,:) = zeros(1,2001);
        stims(count,:) = powers(i)*template/10;
        if ~isempty(spike_times{i}{j})
            responses(count,ceil(spike_times{i}{j})) = 1;
        end
        count = count + 1;
    end
    
end

%%
clc
betahat = fit_lifglm_oneloc(responses,stims,params);

