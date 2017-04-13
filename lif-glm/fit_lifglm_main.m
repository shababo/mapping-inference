% params.g = .01;
% params.current_template = template;
clear responses
clear stims
% loc_ind = 1;
clear stims_ind
count = 1;
spike_count_means = zeros(9,5);
spike_time_means = zeros(9,5);

for k = 1:10%1:length(cells(1).spike_data)
%     if k ~= 1
%         k_ind = k -1;
%     else
%         k_ind = k;
%     end
    spike_times = cells(2).spike_data(k).spike_times;
    powers = cells(1).spike_data(k).powers;
    distances(k_ind) = norm(cells(1).spike_data(k).location(1:2));
    for i = 1:length(spike_times)
        spike_count_means(k_ind,i) = length([spike_times{i}{:}])/length(spike_times{i});
        spike_time_means(k_ind,i) = mean([spike_times{i}{:}]);
        for j = 1:length(spike_times{i})
            responses(count,:) = zeros(1,2001);
            stims(count,:) = powers(i)*template;
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
g = [.009 .01 .03 .05 .07];
devs = zeros(size(g));
for i = 1:length(g)
    g(i)
    params.g = g(i);
    [glm_out(i).beta,fmin_out(i).beta,glm_out(i).stats_conv,fmin_out(i).ll] = ...
        fit_lifglm_oneloc(responses,stims,stims_ind,params);
    devs(i) = glm_out(i).stats_conv.dev;
end



%%
clc
scale = 1;
[min_dev, g_choice] = min(devs);%length(g);
params_sim.V_th = glm_out(g_choice).beta(1)*scale;
params_sim.V_reset = glm_out(g_choice).beta(2)*scale;
num_sim_trials = 100;
params_sim.g = g(g_choice);
funcs.invlink = @invlink_test;%@(resp) log(1 + exp(resp));%@(x) exp(x);
spike_count_means_sim = zeros(9,5);
spike_time_means_sim = zeros(9,5);
for k = [1:9]%1:length(cells(1).spike_data)
%     if k ~= 1
%         k_ind = k -1;
%     else
%         k_ind = k;
%     end

    powers = cells(1).spike_data(k).powers;
    params_sim.gain = glm_out(g_choice).beta(k+2)*scale;
    for j = 1:length(powers)
        spike_times = [];
        for i = 1:num_sim_trials
            
            [V_vect, spikes] = lif_glm_sim(stims((j-1)*5+1,:),params_sim,funcs);
            spike_times = [spike_times find(spikes,1,'first')];
        end
        spike_count_means_sim(k,j) = length(spike_times)/num_sim_trials;
        spike_time_means_sim(k,j) = mean(spike_times);
    end    
end
%%
figure
for i = 1:9
    
    subplot(9,2,(i-1)*2+1)
    plot(powers,spike_count_means_sim(i,:))
    hold on
    plot(powers,spike_count_means(i,:),'--')
    if i == 1
        legend({'sim','data'})
    end
    xlim([0 175])
    subplot(9,2,i*2)
    plot(powers,spike_time_means_sim(i,:))
    hold on
    plot(powers,spike_time_means(i,:),'--')
    xlim([0 175])
end
% params_sim.gain = glm_out(end).beta(3)*scale;
% [V_vect, spikes] = lif_glm_sim(stims(21-10,:),params_sim,funcs);
% figure; plot(V_vect)
% ylim([0-.1*params_sim.V_th params_sim.V_th+.1*params_sim.V_th])

