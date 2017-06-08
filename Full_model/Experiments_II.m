cd('C:/Users/Shizhe/Documents/GitHub/mapping-inference/Full_model/')
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Select cells:
% connected_cells = gamma_ini>0.1;
% %% Pre-calculation
% % Note: use the standard template for inference when testing robustness
% Z_new = Z_dense(connected_cells,:);
% n_cell_temp=sum(connected_cells);
%
% %%
% cell_params.locations = local_locations(connected_cells,:);
% cell_params.shape_gain = ones(n_cell_temp,1);
% shape_template = struct();
% shape_template.shape= l23_average_shape;
% [pi_dense_local, inner_normalized_products] = get_weights_v2(cell_params, shape_template,Z_new);
% %%
% cell_params.locations = local_locations;
% cell_params.shape_gain = local_shape_gain';
% shape_template = l23_cells_for_sim;
% [pi_dense_all, ~] = get_weights_v2(cell_params, shape_template,Z_new);
%
%
% %%
% num_sources=1;
% N2=50; % Number of batches
% grid_index = 1:size(pi_dense_all,2);
%
% % Parameters
% sqrt_transform = false; % whether to use squared transformation
% % Parameters for the working model
% num_threshold=10; % number of bins to use
% mark = 0; % 0: amplitude; 1: latency.
% obj_function = @joint_sparsity_weight_entropy; %
% num_trials_first =max(200, ceil(n_cell_local/num_sources)); % Number of trials in the first batches
% k_minimum = 0.001; % minimum stimulus intensity to consider
% %%
% %---------------------------------------------------------------------%
% % Design stage
% % initialization
% output= struct([]);
% for j = 1:num_threshold
%     output(j).alpha = .1*ones(n_cell_temp*num_power_level+1,1);
%     output(j).mu = zeros(n_cell_temp*num_power_level+1,1);
%     output(j).s_sq = ones(n_cell_temp*num_power_level+1,1);
%     output(j).threshold = [];
% end
% X_g = zeros(0,n_cell_temp*num_power_level);
% locations_trials = zeros(0,num_sources);
% powers_trials= zeros(0,num_sources);
% Y_g = zeros(0,num_threshold);
% counts_freq = zeros(size(Z_new,1),1);
% %%
% cell_params.V_th = 15*ones([n_cell_local,1]);
% cell_params.V_reset = v_reset_known*ones([n_cell_local,1]);
% cell_params.gamma = local_gamma;
% cell_params.amplitudes = abs(normrnd(3,1,[n_cell_local, 1]));
% cell_params.sigma_across =  abs(normrnd(0.5,0.1,[n_cell_local, 1]));
% cell_params.sigma_within = abs(normrnd(0.5,0.1,[n_cell_local, 1]));
% cell_params.locations = local_locations;
% cell_params.shape_gain = local_shape_gain;
% bg_params.mean = 0;
% bg_params.sigma = 1.5;
% bg_params.firing_rate = 0; %spike/sec
% shape_template= l23_cells_for_sim;
% %%
% for i_batch= 1:N2
%     tic
%     tstart=toc;
%     output_old=output;
%     [locations_this_batch, powers_this_batch,counts_freq] = optimal_design_v2(...
%         i_batch, num_sources,num_peaks,num_trials_first,num_trials_batch, output, Y_g, ...
%         num_power_level,random_prop, counts_freq, pi_dense_local,inner_normalized_products,...
%         grid_index, freq_pen, num_samples);
%
%     locations_trials = [locations_trials; locations_this_batch];
%     powers_trials = [powers_trials; powers_this_batch];
%
%     [mpp_temp, ~, ~] = generate_data_v2(...
%         locations_this_batch,powers_this_batch,pi_dense_all,k_minimum,cell_params, shape_template, ...
%         power_level,I_e_vect, data_params,bg_params,trials_specific_variance);
%
%     if i_batch == 1
%         mpp= mpp_temp;
%     else
%         mpp( ((i_batch-2)*num_trials_batch + num_trials_first) + (1:num_trials_batch)) =mpp_temp;
%     end
% end
%
%
%  %%
% stimuli_size_local=zeros(length(mpp),n_cell_temp);
% for l = 1:length(mpp)
%     for m = 1:size(locations_trials,2)
%         stimuli_size_local(l,:)  = stimuli_size_local(l,:)+( pi_dense_local(:,locations_trials(l,m)).*power_level(powers_trials(l,m)))';
%     end
% end
% n_trial = size(stimuli_size_local,1);
%
% evoked_cell = cell(n_trial,1);
% for i_trial = 1:n_trial
%     evoked_cell_index = 0; % 0: background evnets
%     for i_cell = 1:n_cell_temp
%         k = stimuli_size_local(i_trial, i_cell);
%         if k > k_minimum
%             evoked_cell_index = [evoked_cell_index i_cell];
%         end
%     end
%     evoked_cell{i_trial} = evoked_cell_index;
% end
%
%
%
% stimuli_size_temp = stimuli_size_local;
% %%
%
% cell_params.V_th = 15*ones(sum(connected_cells),1);
% cell_params.V_reset = -1e4*ones(sum(connected_cells),1);
%
% cell_params.gain = zeros(n_cell_temp,1);
% for i_cell = 1 : n_cell_temp
%     %cell_params.gain(i_cell) = l23_cells_for_sim(local_shape_gain(i_cell)).optical_gain;
%     cell_params.gain(i_cell) = mean([l23_cells_for_sim.optical_gain]);
%     cell_params.gain_sd(i_cell)= std([l23_cells_for_sim.optical_gain]);
% end
%
% % The local g:
% cell_params.g = zeros(n_cell_temp,1);
% for i_cell = 1 : n_cell_temp
%     %cell_params.g(i_cell) = l23_cells_for_sim(local_shape_gain(i_cell)).g;
%     cell_params.g(i_cell) =  mean([l23_cells_for_sim.g]);
% end
% %
% [estimated_intensity]=Intensity_v5(stimuli_size_temp,  n_stimuli_grid,n_grid_voltage,...
%     t_vect,t_factor,k_minimum,...
%     cell_params, funcs,...
%     I_stimuli,sd_range,delay_params);
%
%
% expected_all = zeros(n_trial,n_cell_temp);
% for i_trial = 1:n_trial
%     for i_cell = 1:n_cell_temp
%         expected_all(i_trial,i_cell)=sum(estimated_intensity{i_trial, i_cell} );
%     end
% end
%
% %% Reformat the mpp:
% for i = 1:n_trial
%    mpp(i).amp=mpp(i).amplitudes;
%    mpp(i).times=mpp(i).event_times;
%
% end

%% Select the cells with decent gammas at the initial fits
%gamma_threshold = 0;
selected_cells = stimulated_cells;
%gamma_initial >gamma_threshold;
%
stimuli_size_selected = stimuli_size_local(:,selected_cells);
n_cell_selected = sum(selected_cells);

evoked_cell_selected = cell(n_trial,1);
for i_trial = 1:n_trial
    evoked_cell_index = 0; % 0: background evnets
    for i_cell = 1:n_cell_selected
        k = stimuli_size_selected(i_trial, i_cell);
        if k > k_minimum
            evoked_cell_index = [evoked_cell_index i_cell];
        end
    end
    evoked_cell_selected{i_trial} = evoked_cell_index;
end

% 1, Initial fits
%   - Estimate the firing rate given initial delay distribution and lif-glm
%   parameters
%   - Estimate the soft assignments and gammas given the fitted values

V_threshold = -50;
cell_params.V_th = 15*ones(n_cell_selected,1);
cell_params.V_reset = v_reset_known*ones(n_cell_selected,1);
cell_params.gain = ones(n_cell_selected,1)*mean([l23_cells_for_sim.optical_gain]);
cell_params.gain_sd= ones(n_cell_selected,1)*std([l23_cells_for_sim.optical_gain]);
cell_params.g =  ones(n_cell_selected,1)*mean([l23_cells_for_sim.g]);

n_delay_grid = 200;
outputM=false;
delay_params_est.type=1;
delay_params_est.mean=35*ones(n_cell_selected,1);
delay_params_est.std=5*ones(n_cell_selected,1);


%% Iterative updates with filtered cells 

gain_fits=zeros(n_cell_selected,maxit);
gamma_fits=zeros(n_cell_selected,maxit);
mu_fits=zeros(n_cell_selected,maxit);
delay_mean_fits=zeros(n_cell_selected,maxit);
delay_std_fits=zeros(n_cell_selected,maxit);

gain_old = 0.003*ones(n_cell_selected,1);
gamma_old= 0.009*ones(n_cell_selected,1);
mu_old = 2*ones(n_cell_selected,1);
sigma_old = ones(n_cell_selected,1);


sparsity_params.threshold=4;
sparsity_params.eta=0.1;
sparsity_params.sparsity=1;

soft_threshold=0.1;
num_MC_lifglm = 5;

normalized_change_outer=1;
convergence_epsilon_outer=0.01;
num_iter=1;
maxit=100;

% Obtain delay vector
delay_params_est.type=1;
delay_params_est.mean=35*ones(n_cell_selected,1);
delay_params_est.std=10*ones(n_cell_selected,1);

n_delay_grid = 200;

%%
while (normalized_change_outer > convergence_epsilon_outer) & (num_iter < maxit)
    num_iter = num_iter+1;
    
    [Stimuli_grid, Intensity_grid]=Intensity_v8(...
        stimuli_size_selected, mpp,I_stimuli,... % data from exp
        cell_params,... % estimated parameters
        funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
        n_stimuli_grid,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
        V_threshold,stimulus_threshold,first_only);
    
    % Get the event rates and expec ted event counts
    expected_all = zeros(n_trial,n_cell_selected);
    event_rates = cell(n_trial,1);
    for i_trial = 1:n_trial
        if length(mpp(i_trial).times)>0
            event_rates{i_trial}=zeros(length(mpp(i_trial).times),n_cell_selected);
        end
    end
    for i_cell = 1:n_cell_selected
        
        %-----------------------------------------%
        % Convolute the intensity with delay distribution 
        delay_prob = zeros(2*n_delay_grid+1,1);
        if delay_params.std == 0
            delay_prob(n_delay_grid+1)=1;
            min_delay=0;
            max_delay=0;
        else
            delay_prob = normpdf( -(-n_delay_grid:n_delay_grid),...
                delay_params_est.mean(i_cell),delay_params_est.std(i_cell));
            % we approximate the probability with densities
            delay_prob = delay_prob/sum(delay_prob);
            min_delay = 1-1-n_delay_grid;
            max_delay = length(delay_prob)-1-n_delay_grid;
        end
        for i_stimuli = 1:length(Stimuli_grid)
            M_grid_intensity{i_stimuli}=zeros(length(Intensity_grid{1}),1);
            for i_t = 1:length(Intensity_grid{1})
                idx_time = max(i_t+min_delay,1): min(i_t+max_delay,n_grid_time);
                idx_delay = (min(idx_time)-i_t+n_delay_grid+1) : (max(idx_time)-i_t+n_delay_grid+1);
                M_grid_intensity{i_stimuli}(i_t)=delay_prob(idx_delay)*Intensity_grid{i_stimuli}(idx_time);
            end
        end
        %------------------------------------------------%
        for i_trial = 1:n_trial
           %------------------------------------------------%
            % Allowing for noise in the gain estimates 
            k_up = stimuli_size_selected(i_trial, i_cell)*(cell_params.gain(i_cell)+sd_range*cell_params.gain_sd(i_cell));
            k_low = stimuli_size_selected(i_trial, i_cell)*(cell_params.gain(i_cell)-sd_range*cell_params.gain_sd(i_cell));
            intensity_temp = zeros([length(t_vect) 1]);
            index_seq = (k_low<Stimuli_grid &  k_up>Stimuli_grid);
            if sum(index_seq)>0
                for i_grid = 1:length(Stimuli_grid)
                    if index_seq(i_grid)>0
                        intensity_temp= intensity_temp+M_grid_intensity{i_grid};
                    end
                end
                intensity_temp=intensity_temp/sum(index_seq);
                expected_all(i_trial,i_cell)=sum(intensity_temp);
                if length(mpp(i_trial).times)>0
                    event_rates{i_trial}(:,i_cell)=intensity_temp( max(1,round(mpp(i_trial).times)) );
                end
            end
            %------------------------------------------------%
        end
       % fprintf('%d\n',i_cell);
    end
    %------------------------------------------------%
    % Updating the gammas and the soft assignments
    [gamma_path mu_path sigma_path total_time soft_assignments bg_rate]= ...
        EM_fullmodel_v3(mpp(1:n_trial), ...
        event_rates,...
        evoked_cell_selected,expected_all, ...
        n_cell_local, gamma_old, mu_old, sigma_old, ...
        convergence_epsilon,f_background, mean_background, sigma_background, ...
        maxit,t_vect,use_size,background_update,sparsity_params);
    %-----------------------------------------------%
    
    
    fprintf('Changes %d\n', normalized_change_outer);
    fprintf('Sum of gamma %d\n', sum( gamma_path(:,end)));
    fprintf('Sum of gain %d\n', sum( gain_old));
    
    % lif-glm updates:
    % Use Monte-Carlo method to update lif-glm parameters based on soft
    % assignments
    % should be turned into a function
    
    % Reformat the soft assignments for each cell
    
    soft_assignments_by_cell = cell(n_cell_selected,1);
    % record the trial, spike time, and soft assignments
    for i_trial = 1:n_trial
        n_event = length(mpp(i_trial).event_times);
        cell_list = evoked_cell_selected{i_trial};
        if n_event >0
            if length(cell_list)>1 % at least one cell beside the background
                for i_cell = 2:length(cell_list)
                    if length(soft_assignments_by_cell{cell_list(i_cell)})==0
                        soft_assignments_by_cell{cell_list(i_cell)}=[i_trial*ones(1,n_event); mpp(i_trial).event_times;...
                            soft_assignments{i_trial}(:,i_cell)'];
                    else
                        soft_assignments_by_cell{cell_list(i_cell)}=[soft_assignments_by_cell{cell_list(i_cell)}...
                            [i_trial*ones(1,n_event); mpp(i_trial).event_times;soft_assignments{i_trial}(:,i_cell)']];
                    end
                end
            end
        end
    end
    
    gains_sample = zeros(n_cell_selected,num_MC_lifglm);
    gains_sd_sample= zeros(n_cell_selected,num_MC_lifglm);
    delay_params_sample.mean = zeros(n_cell_selected,num_MC_lifglm);
    delay_params_sample.std= zeros(n_cell_selected,num_MC_lifglm);
    cell_data = cell(n_cell_selected,1);
    
    % Calculate the crude grid with margins of errors 
    for i_MC = 1:num_MC_lifglm
        %-------------------------------------------------%
        % Draw the hard assignments 
        % Draw one event for each cell
%         t1=toc;
        for i_cell = 1:n_cell_selected
            soft_temp =soft_assignments_by_cell{i_cell};
             cell_data{i_cell}=struct();
                cell_data{i_cell}.responses = [];
                cell_data{i_cell}.stims = [];
            if length(soft_temp)<1
               % Do nothing..
            else
                trial_list=  unique(soft_temp(1,:));
                for i_trial = 1:length(trial_list)
                    trial_idx= soft_temp(1,:)==trial_list(i_trial);
                    soft_this_trial = soft_temp(:,trial_idx);
                    prob_sum = sum(soft_this_trial(3,:));
                    response_this_trial = zeros(1,length(I_stimuli));
                    if prob_sum > soft_threshold
                        if prob_sum > 1
                            prob_sample = soft_this_trial(3,:)/prob_sum;
                        else
                            prob_sample=soft_this_trial(3,:);
                        end
                        prob_sample = [1- sum(prob_sample) prob_sample];
                        r_temp = rand(1);
                        i_event = min(find(r_temp<cumsum(prob_sample)));
                        if i_event > 1
                            response_this_trial(max(1,round(soft_this_trial(2,i_event-1))))=1;
                        end
                    end
                    if sum(response_this_trial)>0
                        cell_data{i_cell}.responses = [cell_data{i_cell}.responses; response_this_trial];
                        cell_data{i_cell}.stims =  [cell_data{i_cell}.stims; ...
                            I_stimuli*stimuli_size_selected(trial_list(i_trial), i_cell)];
                    end
                end
            end
        end
%         t2=toc;
%         timevect(1)=t2-t1;
        %----------------------------------------------------------------%
        
        %----------------------------------------------------------------%
        % Update the LIFGLM parameters and the delay distributions 
        delays=[];
                
        lif_glm_gains= zeros(n_cell_selected,1);
        delay_params_temp.mean = zeros(n_cell_selected,1);
        delay_params_temp.std = zeros(n_cell_selected,1);
        for i_cell = 1:n_cell_selected
            N_cell = length(cell_data{i_cell}.responses);
            in_params.g = cell_params.g(i_cell);
            if (N_cell/length(I_stimuli))>0
                responses=cell_data{i_cell}.responses;
                stims=cell_data{i_cell}.stims;
                
                n_trial_temp = size(responses,1);
                responses_reg=responses;responses_reg(:,:)=0;
                for i_trial = 1:n_trial_temp
%                      t3=toc;
                    k_temp = max(stims(i_trial,:))/max(I_stimuli);
                    k_up =k_temp*(cell_params.gain(i_cell)+sd_range*cell_params.gain_sd(i_cell));
                    k_low =  k_temp*(cell_params.gain(i_cell)-sd_range*cell_params.gain_sd(i_cell));
                    intensity_temp = zeros([length(t_vect) 1]);
                    index_seq = (k_low<Stimuli_grid &  k_up>Stimuli_grid);
                    if sum(index_seq)>0
                        for i_grid = 1:length(Stimuli_grid)
                            if index_seq(i_grid)>0
                                intensity_temp= intensity_temp+M_grid_intensity{i_grid};
                            end
                        end
                        intensity_temp=intensity_temp/sum(index_seq);
                    end
%                     t4=toc;
                   [~, idx_min]=min( abs(k_temp-Stimuli_grid));
                   %intensity_temp=M_grid_intensity_error{idx_min};
                    % Find the MLE delay given the intensity and the delay
                    % distribution 
                    spikes=find(responses(i_trial,:));
                    if length(spikes)>0
                        spike_first = spikes(1);
                        delay_seq =spike_first-(1:length(I_stimuli));
                        delay_prob = normpdf(delay_seq,...
                            delay_params_est.mean(i_cell),delay_params_est.std(i_cell));
                        prod_prob = delay_prob.*intensity_temp';
                        [~, spike_max]= max(prod_prob);
                        responses_reg(i_trial,spike_max)=1;
                        delays =[delays spike_first-spike_max];
                    end
                    % Update the delay distribution 
%                     delay_params_temp.mean(i_cell)=mean(delays);
%                     if std(delays)==0
%                        % not updating the standard deviations  
%                     else
%                         delay_params_temp.std(i_cell)=delay_params_temp.std(i_cell);
%                     end
%                     t5=toc;
%                     timevect(2)=timevect(2)+t5-t4;
                end
                    % Fit the LIF-GLM using the adjusted spikes
                    %-------------------------------------%
                    %lif_glm_gains(i_cell)=stats_conv.beta(2);
                 [stats_conv] = fit_lifglm_v3(responses_reg, stims,in_params,v_reset_known,first_only);
%                     t6=toc;
%                     timevect(3)=timevect(3)+t6-t5;
                    lif_glm_gains(i_cell)=stats_conv.beta;
            else
                lif_glm_gains(i_cell)=cell_params.gain(i_cell);
%                 delay_params_temp.mean(i_cell)=delay_params_est.mean(i_cell);
%                 delay_params_temp.std(i_cell)=delay_params_est.std(i_cell);
            end
        end
        %delay_params_sample.mean(:,i_MC) = delay_params_temp.mean;
        %delay_params_sample.std(:,i_MC)= delay_params_temp.std;
        delay_params_sample.mean(:,i_MC) = mean(delays)*ones(n_cell_selected,1);
        delay_params_sample.std(:,i_MC)= std(delays)*ones(n_cell_selected,1);
        
        gains_sample(:,i_MC)=lif_glm_gains;
        fprintf('%d MC sample completed\n',i_MC);
    end
    
    % Evaluate the updates:
    gamma_current= gamma_path(:,end);
    sigma_current = sigma_path(:,end);
    mu_current = mu_path(:,end);
    gain_current = median(gains_sample,2);
   
    for i_cell = 1:n_cell_selected
        gain_sd_current(i_cell) = std(gains_sd_sample(i_cell,:));
        if gain_sd_current(i_cell)==0
            gain_sd_current(i_cell)=cell_params.gain_sd(i_cell);
        end
    end
    
    delay_params_est.mean = mean(mean(delay_params_sample.mean,2))*ones(n_cell_selected,1);
    %mean(delay_params_sample.mean,2);
    %
    delay_params_est.std= mean(mean(delay_params_sample.std,2))*ones(n_cell_selected,1);
    %mean(delay_params_sample.std,2);
%     delay_params_est.mean = 35*ones(n_cell_selected,1);
%     delay_params_est.std=10*ones(n_cell_selected,1);
   
    normalized_change_outer = norm(gamma_current - gamma_old)/(norm(gamma_old)+1) + norm(mu_current - mu_old)/(norm(mu_old)+1)+...
        norm(sigma_current - sigma_old)/(norm(sigma_old)+1)+norm(gain_current-gain_old)/(norm(gain_old)+1);
    
    gamma_old =  gamma_current;
    sigma_old = sigma_current;
    mu_old = mu_current;
    gain_old = gain_current;
    
    
    % Update the intensities
    cell_params.gain = gain_current;
    cell_params.gain_sd = gain_sd_current;
    
    gain_fits(:,num_iter)=gain_current;
    mu_fits(:,num_iter)=mu_current;
    gamma_fits(:,num_iter)=gamma_current;
    delay_mean_fits(:,num_iter)=delay_params_est.mean;
    delay_std_fits(:,num_iter)=delay_params_est.std;
    
    
figure(num_iter)

         plot(local_gamma(stimulated_cells)+normrnd(0,0.01,[sum(stimulated_cells) 1]),...
             gamma_initial(stimulated_cells),'.','MarkerSize',20,'col','b')
         hold on;
     xlim([-0.1,1.1]);
ylim([-0.1,1.1]);
xlabel('True gamma')
ylabel('Est. gamma')
gamma_current = gamma_fits(:,num_iter);
gamma_all = zeros(n_cell_local,1);
gamma_all(selected_cells)=gamma_current;
plot(local_gamma(selected_cells)+normrnd(0,0.01,[sum(selected_cells) 1]),...
    gamma_all(selected_cells),'.','MarkerSize',20,'col','r')
xlim([-0.1,1.1]);
ylim([-0.1,1.1]);

line([0 1],[0 1]);
xlabel('True gamma')
ylabel('Est. gamma')
hold off;
end
%%

figure(1)

    gamma_initial =zeros(n_cell_local,1);
    gamma_initial(stimulated_cells)=gamma_ini;
         plot(local_gamma(stimulated_cells)+normrnd(0,0.01,[sum(stimulated_cells) 1]),...
             gamma_initial(stimulated_cells),'.','MarkerSize',20,'col','b')
         hold on;
     xlim([-0.1,1.1]);
ylim([-0.1,1.1]);
xlabel('True gamma')
ylabel('Est. gamma')

gamma_current = gamma_fits(:,num_iter-1);
gamma_all = zeros(n_cell_local,1);
gamma_all(selected_cells)=gamma_current;
plot(local_gamma(selected_cells)+normrnd(0,0.01,[sum(selected_cells) 1]),...
    gamma_all(selected_cells),'.','MarkerSize',20,'col','r')
xlim([-0.1,1.1]);
ylim([-0.1,1.1]);

line([0 1],[0 1]);
xlabel('True gamma')
ylabel('Est. gamma')
hold off;
%%
gain_all = zeros(n_cell_local,1);
gain_all(selected_cells)=gain_current;


plot(local_gain(local_gamma>0)+normrnd(0,0.001,[sum(local_gamma>0) 1]),gain_all(local_gamma>0),'.','MarkerSize',20)
xlim([0,0.06]);
ylim([0,0.06]);
line([0 1],[0 1]);
xlabel('True gain')
ylabel('Est. gain')

