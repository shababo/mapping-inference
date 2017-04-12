params.g = .05;
% params.current_template = template;
clear responses
clear stims
% loc_ind = 1;

count = 1;
for k = 1:length(cells(1).spike_data)
    spike_times = cells(1).spike_data(k).spike_times;
    powers = cells(1).spike_data(k).powers;
    distances(k) = norm(cells(1).spike_data(k).location(1:2));
    for i = 1:length(spike_times)
        
        for j = 1:length(spike_times{i})
            responses(count,:) = zeros(1,2001);
            stims(count,:) = powers(i)*template/1000;
            stims_ind(count) = k;
            if ~isempty(spike_times{i}{j})
                responses(count,ceil(spike_times{i}{j})) = 1;
            end
            count = count + 1;
        end
    end    
end


%%
clc
betahat = fit_lifglm_oneloc(responses,stims,stims_ind,params);
%%
figure; scatter(distances,betahat(3:end))
%%


