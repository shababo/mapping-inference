params.g = .01;
% params.current_template = template;
clear responses
clear stims
% loc_ind = 1;
clear stims_ind
count = 1;
for k = [1 3:10]%1:length(cells(1).spike_data)
    if k ~= 1
        k_ind = k -1;
    else
        k_ind = k;
    end
    spike_times = cells(1).spike_data(k).spike_times;
    powers = cells(1).spike_data(k).powers;
    distances(k_ind) = norm(cells(1).spike_data(k).location(1:2));
    for i = 1:length(spike_times)
        
        for j = 1:length(spike_times{i})
            responses(count,:) = zeros(1,2001);
            stims(count,:) = powers(i)*template/1000;
            stims_ind(count) = k_ind;
            if ~isempty(spike_times{i}{j})
                responses(count,ceil(spike_times{i}{j})) = 1;
            end
            count = count + 1;
        end
    end    
end


%%
clc
clear glm_out
clear fmin_out
g = [1e-1 1e-2];
for i = 1:length(g)
    params.g = g(i);
    [glm_out(i).beta,fmin_out(i).beta,glm_out(i).stats_conv,fmin_out(i).ll] = ...
        fit_lifglm_oneloc(responses,stims,stims_ind,params);
end
%%
figure; scatter(distances,betahat_fmin(3:end)); hold on;
scatter(distances,betahat_glm(4:end)); 
%%
clc
scale = 100;
params_sim.V_th = glm_out(end).beta(1)*scale;
params_sim.V_reset = glm_out(end).beta(2)*scale;
params_sim.gain = glm_out(end).beta(3)*scale;
params_sim.g = g(end);
funcs.invlink = @invlink_test;%@(resp) log(1 + exp(resp));%@(x) exp(x);

[V_vect, spikes] = lif_glm_sim(stims(21-0,:),params_sim,funcs);
figure; plot(V_vect)
% ylim([0-.1*params_sim.V_th params_sim.V_th+.1*params_sim.V_th])

