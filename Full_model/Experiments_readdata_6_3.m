% Loading functions and Data generation
% cd('C:/Users/Shizhe/Documents/GitHub/mapping-inference/Full_model/')
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Paramters in the simulations
trials_specific_variance= 0;
load('./Environments/l23_template_cell.mat');
%l23_average_shape
temp=l23_average_shape;
temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;

load('./Environments/l23_cells_for_sim.mat');
% We only use the main gain in this analysis 
mean_gain = mean([l23_cells_for_sim.optical_gain]);
%% load real data
load('./Environments/6_3_s2c3_mpp_and_stim_data.mat')
%cell_locs, target_locs,

%% Setting experiment details based on the data
local_locations= cell_locs;
Z_dense =target_locs;
powers_trials = stim_pow;
locations_trials=target_inds;
%% Set seed for reproducibility
rng(12242,'twister');
n_trial = size(locations_trials,1);
% n_trial = 8000;
n_cell_local = size(local_locations,1);

%% Calculate the background rate 
mpp_copy=mpp;
histogram([mpp_copy.times]);
% Why are there a lot of events around 500?

trial_counts = 0;
event_counts = 0;
for i_trial = 1:n_trial
    %if stim_pow(i_trial)==50
    trial_counts = trial_counts+1;
    event_counts = event_counts +sum(mpp_copy(i_trial).times>0 & mpp_copy(i_trial).times<80);
    %end
end
background_rate = event_counts/trial_counts/80;
%%
cell_params.locations =  local_locations;
cell_params.shape_gain = ones(n_cell_local,1);
shape_template = struct();
shape_template.shape= l23_average_shape;
[pi_dense_local, inner_normalized_products] = get_weights_v2(cell_params, shape_template,Z_dense);
%% Loading the current template using new template
load('./Environments/chrome-template-3ms.mat');
downsamp=1;
max_time=300;

power_level = unique(powers_trials);
num_power_level=length(power_level);
current_template=template(1:downsamp:max_time);

I_stimuli =current_template;
evoked_params.stim_start = 1;
evoked_params.stim_end = length(I_stimuli);
data_params.T = length(I_stimuli); % total time at 20k Hz
data_params.dt = 1; % 1/20 ms

maxit=100;
convergence_epsilon=1e-3;
mean_background = 2;
sigma_background  = 1;
%gamma_threshold=0.1;
k_minimum = 0.1; % minimum stimulus intensity to consider

%%
v_reset_known=-4e3;
cell_params.V_th = 15*ones([n_cell_local,1]);
cell_params.V_reset = v_reset_known*ones([n_cell_local,1]);

%%
    max_time = 300;
min_time=60;
for i_trial = 1:n_trial
    if mpp(i_trial).num_events >0
        range_idx = mpp(i_trial).times<max_time & mpp(i_trial).times>min_time ;

        mpp(i_trial).num_events = sum(range_idx);
        mpp(i_trial).times = mpp(i_trial).times(range_idx);
        mpp(i_trial).amp =mpp(i_trial).amp(range_idx);
    end

end
%
%     %% count the number of events:
%    sum([mpp.assignments]==0)
 %   240*n_trial*background_rate
%
%
%     length([mpp.assignments])
%  length([mpp_copy(1:n_trial).times])
%

% save('mpp_sim_II.mat','mpp');
%%
%------------------------------------%
% Estimating the marginal firing rate
sigma_unknown=1;
T=length(I_stimuli); % total time at 20k Hz
dt=1;
t_vect= dt:dt:T;
sd_range=1.5;
n_stimuli_grid=20;
n_grid_voltage=400;
t_factor=1;
gap_stimuli=0.5;
first_only=true;
stimulus_threshold=0.1;
funcs.invlink = @invlink_test;%@(resp) log(1 + exp(resp));%@(x) exp(x);
n_grid_time = length(I_stimuli);
%%

stimuli_size_local=zeros(n_trial,n_cell_local);
for l = 1:n_trial
    for m = 1:size(locations_trials,2)
        if isnan(locations_trials(l,m))
        else
            stimuli_size_local(l,:) = stimuli_size_local(l,:)+( pi_dense_local(:,locations_trials(l,m)).*powers_trials(l))';
        end
    end
end
n_trial = size(stimuli_size_local,1);
%% Throw away cells that are not stimulated enough:
stim_threshold = 20;
stimulated_cells = sum(stimuli_size_local>stim_threshold )>5;
%%
stimuli_size_stimulated = stimuli_size_local(:,stimulated_cells);
n_cell_stimulated = sum(stimulated_cells);

evoked_cell_stimulated = cell(n_trial,1);
for i_trial = 1:n_trial
    evoked_cell_index = 0; % 0: background evnets
    for i_cell = 1:n_cell_stimulated
        k = stimuli_size_stimulated(i_trial, i_cell);
        if k > k_minimum
            evoked_cell_index = [evoked_cell_index i_cell];
        end
    end
    evoked_cell_stimulated{i_trial} = evoked_cell_index;
end


%% Analysis:
% 1, Initial fits
%   - Estimate the firing rate given initial delay distribution and lif-glm
%   parameters
%   - Estimate the soft assignments and gammas given the fitted values

V_threshold = -50;
cell_params.V_th = 15*ones(sum(stimulated_cells),1);
cell_params.V_reset = v_reset_known*ones(sum(stimulated_cells),1);
cell_params.gain = ones(n_cell_stimulated,1)*mean([l23_cells_for_sim.optical_gain]);
cell_params.gain_sd= ones(n_cell_stimulated,1)*std([l23_cells_for_sim.optical_gain]);
cell_params.g =  ones(n_cell_stimulated,1)*mean([l23_cells_for_sim.g]);

n_delay_grid = 200;
outputM=false;
delay_params_est.type=1;
delay_params_est.mean=35*ones(n_cell_stimulated,1);
delay_params_est.std=10*ones(n_cell_stimulated,1);


[Stimuli_grid, Intensity_grid]=Intensity_v8(...
    stimuli_size_stimulated, mpp,I_stimuli,... % data from exp
    cell_params,... % estimated parameters
    funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
    n_stimuli_grid,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
    V_threshold,stimulus_threshold,first_only);

expected_all = zeros(n_trial,n_cell_stimulated);
event_rates = cell(n_trial,1);
for i_trial = 1:n_trial
    if length(mpp(i_trial).times)>0
        event_rates{i_trial}=zeros(length(mpp(i_trial).times),n_cell_stimulated);
    end
end
for i_cell = 1:n_cell_stimulated
    
    %-----------------------------------------%
    % Convolute the intensity with delay distribution
    delay_prob = zeros(2*n_delay_grid+1,1);
    if delay_params_est.std(i_cell) == 0
        delay_prob(n_delay_grid+1)=1;
        min_delay=0;
        max_delay=0;
    else
        delay_prob = normpdf( -(-n_delay_grid:n_delay_grid),delay_params_est.mean(i_cell),delay_params_est.std(i_cell));
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
        k_up = stimuli_size_stimulated(i_trial, i_cell)*(cell_params.gain(i_cell)+sd_range*cell_params.gain_sd(i_cell));
        k_low = stimuli_size_stimulated(i_trial, i_cell)*(cell_params.gain(i_cell)-sd_range*cell_params.gain_sd(i_cell));
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
    fprintf('%d\n',i_cell);
end
%%
background_update=0;
f_background = background_rate;

mean_background = 0; %not actually used
sigma_background  =1;%not actually used

gain_old = 0.003*ones(n_cell_stimulated,1);
gamma_old= 0.009*ones(n_cell_stimulated,1);

mu_old = 2*ones(n_cell_stimulated,1);
sigma_old = ones(n_cell_stimulated,1);

use_size =0;
convergence_epsilon = 0.01;
maxit = 100;

sparsity_params.threshold=4;
sparsity_params.eta=0.1;
sparsity_params.sparsity=0;

[gamma_path mu_path sigma_path total_time soft_assignments bg_rate]= ...
    EM_fullmodel_v3(mpp(1:n_trial), ...
    event_rates,...
    evoked_cell_stimulated,expected_all, ...
    n_cell_local, gamma_old, mu_old, sigma_old, ...
    convergence_epsilon,f_background, mean_background, sigma_background, ...
    maxit,t_vect,use_size,background_update,sparsity_params);

gamma_ini= gamma_path(:,end);
sigma_ini = sigma_path(:,end);
mu_ini = mu_path(:,end);

%% Visualize the initial fits
gamma_initial =zeros(n_cell_local,1);
gamma_initial(stimulated_cells)=gamma_ini;
save('initial_fits_6_3.mat','gamma_initial');


%% Select the cells with decent gammas at the initial fits
gamma_threshold = 0.05;
selected_cells = gamma_initial >gamma_threshold;
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
sparsity_params.sparsity=0;

soft_threshold=0.1;
num_MC_lifglm = 5;

normalized_change_outer=1;
convergence_epsilon_outer=0.01;
num_iter=1;
maxit=100;

% Initialize the delay distribution:
delay_params_est.type=1;
delay_params_est.mean=35*ones(n_cell_selected,1);
delay_params_est.std=15*ones(n_cell_selected,1);


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
        if delay_params_est.std == 0
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
        n_event = length(mpp(i_trial).times);
        cell_list = evoked_cell_selected{i_trial};
        if n_event >0
            if length(cell_list)>1 % at least one cell beside the background
                for i_cell = 2:length(cell_list)
                    if length(soft_assignments_by_cell{cell_list(i_cell)})==0
                        soft_assignments_by_cell{cell_list(i_cell)}=[i_trial*ones(1,n_event); mpp(i_trial).times;...
                            soft_assignments{i_trial}(:,i_cell)'];
                    else
                        soft_assignments_by_cell{cell_list(i_cell)}=[soft_assignments_by_cell{cell_list(i_cell)}...
                            [i_trial*ones(1,n_event); mpp(i_trial).times;soft_assignments{i_trial}(:,i_cell)']];
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
    
   
end
%%

gamma_all = zeros(n_cell_local,1);
gamma_all(selected_cells)=gamma_current;
gain_all = zeros(n_cell_local,1);
gain_all(selected_cells)=gain_current;

save('final_fits_6_3.mat','gamma_all','gain_all');

