% params.g = .01;
% params.current_template = template;
clear responses
clear stims
% loc_ind = 1;
clear stims_ind
count = 1;
K = 10;
spike_count_means2 = zeros(K,length(powers));
spike_time_means = zeros(K,length(powers));
cell_choice = 1;
downsamp = 1;
for k = 1:K%1:length(cells(1).spike_data)
%     if k ~= 1
%         k = k -1;
%     else
%         k = k;
%     end
    spike_times = cells(cell_choice).spike_data(k).spike_times;
    powers = cells(cell_choice).spike_data(k).powers;
    distances(k) = norm(cells(cell_choice).spike_data(k).location);
    for i = 1:length(spike_times)
        spike_count_means2(k,i) = length([spike_times{i}{:}])/length(spike_times{i});
        spike_time_means(k,i) = mean(floor([spike_times{i}{:}]/downsamp));
        for j = 1:length(spike_times{i})
            responses(count,:) = zeros(1,length(template(1:downsamp:end)));
            stims(count,:) = powers(i)*template(1:downsamp:end);
            stims_ind(count) = k;
            if ~isempty(spike_times{i}{j})
                responses(count,floor(spike_times{i}{j}/downsamp)) = 1;
            end
            count = count + 1;
        end
    end    
end


%%
clc
clear glm_out
clear fmin_out
% g = [.009 .01 .03 .05 .07];
g = .05*downsamp;
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
num_sim_trials = 50;
params_sim.g = g(g_choice);
funcs.invlink = @invlink_test;%@(resp) log(1 + exp(resp));%@(x) exp(x);
spike_count_means_glmfit_sim = zeros(9,5);
spike_time_means_glmfit_sim = zeros(9,5);
for k = [1:10]%1:length(cells(1).spike_data)
%     if k ~= 1
%         k = k -1;
%     else
%         k = k;
%     end
    k
    powers = cells(1).spike_data(k).powers;
    params_sim.gain = glm_out(g_choice).beta(k+2)*scale;
    for j = 1:length(powers)
        spike_times = [];
        for i = 1:num_sim_trials
            
            [V_vect, spikes] = lif_glm_sim(stims((j-1)*5+1,:),params_sim,funcs);
            spike_times = [spike_times find(spikes,1,'first')];
        end
        spike_count_means_glmfit_sim(k,j) = length(spike_times)/num_sim_trials;
        spike_time_means_glmfit_sim(k,j) = mean(spike_times);
        
    end    
end

% [min_dev, g_choice] = min(devs);%length(g);
% params_sim.V_th = fmin_out(g_choice).beta(1)*scale;
% params_sim.V_reset = fmin_out(g_choice).beta(2)*scale;
% num_sim_trials = 100;
% params_sim.g = g(g_choice);
% funcs.invlink = @invlink_test;%@(resp) log(1 + exp(resp));%@(x) exp(x);
% spike_count_means_fmin_sim = zeros(9,5);
% spike_time_means_fmin_sim = zeros(9,5);
% for k = [1:9]%1:length(cells(1).spike_data)
% %     if k ~= 1
% %         k = k -1;
% %     else
% %         k = k;
% %     end
%     k
%     powers = cells(1).spike_data(k).powers;
%     params_sim.gain = fmin_out(g_choice).beta(k+2)*scale;
%     for j = 1:length(powers)
%         spike_times = [];
%         for i = 1:num_sim_trials
%             
%             [V_vect, spikes] = lif_glm_sim(stims((j-1)*5+1,:),params_sim,funcs);
%             spike_times = [spike_times find(spikes,1,'first')];
%         end
%         spike_count_means_fmin_sim(k,j) = length(spike_times)/num_sim_trials;
%         spike_time_means_fmin_sim(k,j) = mean(spike_times);
%         
%     end    
% end
%%
figure
for i = 1:9
    
    subplot(9,2,(i-1)*2+1)
    plot(powers,spike_count_means(i,:))
    hold on
%     plot(powers,spike_count_means_fmin_sim(i,:),'--')
    hold on
    plot(powers,spike_count_means_glmfit_sim(i,:),'--')
    if i == 1
        legend({'data','fmin','glmfit'})
    end
    xlim([0 175])
    subplot(9,2,i*2)
    plot(powers,spike_time_means(i,:))
    hold on
%     plot(powers,spike_time_means_fmin_sim(i,:),'--')
    hold on
    plot(powers,spike_time_means_glmfit_sim(i,:),'--')
    xlim([0 175])
    ylim([0 25]*20/downsamp)
%     subplot(9,3,i*3)
%     plot(powers,spike_count_means_fmin_sim(i,:))
%     hold on
%     plot(powers,spike_time_means_fmin_sim(i,:),'--')
%     xlim([0 175])
end
% params_sim.gain = glm_out(end).beta(3)*scale;
% [V_vect, spikes] = lif_glm_sim(stims(21-10,:),params_sim,funcs);
% figure; plot(V_vect)
% ylim([0-.1*params_sim.V_th params_sim.V_th+.1*params_sim.V_th])

