%% Testing the effect of delay on lif-glm & firing rate

addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%%
rng(12242,'twister');
num_sample = 500;
load('./Environments/l23_cells_for_sim.mat');
num_types_cell = length(l23_cells_for_sim);
% normalized the cell shapes
for i = 1:num_types_cell
    temp=l23_cells_for_sim(i).shape;
    temp_max = max(max(max(temp)));
    l23_cells_for_sim(i).shape = temp/temp_max;
end
%
i_template = 1;

%
v_reset_known=-4000;
params_sim.V_th= 15;
params_sim.V_reset = v_reset_known;
params_sim.g =  l23_cells_for_sim(i_template).g;
params_sim.gain =  l23_cells_for_sim(i_template).optical_gain;
funcs.invlink = @invlink_test;%@(resp) log(1 + exp(resp));%@(x) exp(x);
stoc_params.mu=0;
stoc_params.sigma = 1e-8;

%stimuli_seq = num_sample*(rand([num_sample 1])+0.5);
stimuli_seq = 50*ones(num_sample,1);

%
load('./Environments/chrome-template-3ms.mat');
downsamp=1;
time_max= 300;
current_template=template(1:downsamp:time_max);
I_e_vect=current_template;
%

% delay_params.mean=30;
% delay_params.std=5;

%delay = normrnd(delay_params.mean,delay_params.std);
%delay = 40;


n_grid_voltage =1000;
stimulus_threshold=-50;
dt=1;
n_grid_time = length(I_e_vect);
I_stimuli= I_e_vect;
t_factor=1;
k_minimum=0.001;

sd_range=1;

%delay_params.mean=round(delay);
%delay_params.std=0.01;
n_delay_grid = 200;

stimulus_threshold=0.1;
n_stimuli_grid=0;
gap_stimuli=5;

stimuli_size=stimuli_seq;

outputM=true;

num_MC=100;
first_only=true;
V_threshold = -50;
%%
for i_outer = 1:3
    responses = zeros(num_sample, length(I_e_vect));
    stims = zeros(num_sample, length(I_e_vect));
    
    cell_params=params_sim;
    cell_params.gain_sd=0.01;
    
    if i_outer==1
        delay_params.type=1;
        delay_params.mean=0;
        delay_params.std=0.1;
    elseif i_outer==2
        delay_params.type=0;
        delay_params.shape=30*30/25;
        delay_params.scale=25/30;
    else
        delay_params.type=1;
        delay_params.mean=35;
        delay_params.std=5;
    end
    for i_trial = 1:num_sample
        k=stimuli_seq(i_trial);
        stim = I_e_vect*k;
        stims(i_trial,:) = stim;
        [V_vect, spikes]  = lif_glm_sim_v2(stim,params_sim,funcs);
        
        % shift the spikes give the delay
        
        if i_outer==1
            delay=0;
        elseif i_outer==2
            delay = round(gamrnd(delay_params.shape,delay_params.scale,[sum(spikes) 1]));
        elseif i_outer==3
            delay = round(normrnd(delay_params.mean,delay_params.std,[sum(spikes) 1]));
        end
        
        
        spikes_delayed=find(spikes)+round(delay)';
        spikes_delayed=spikes_delayed( spikes_delayed<time_max);
        responses(i_trial,spikes_delayed)=1;
        %sum(responses(i_trial,:))
    end
    
    mpp=struct();
    for i=1:num_sample
        mpp(i).times=find(responses(i,:));
    end
    avg_spikes{i_outer}=mean(responses,1);
    
    
    % draw MC samples of real spikes
    for i_MC = 1:num_MC
        responses_MC=zeros(num_sample, length(I_e_vect));
        
        for i_trial = 1:num_sample
            
            % shift the spikes give the delay
            
            if i_outer==1
                delay=0;
            elseif i_outer==2
                delay = round(gamrnd(delay_params.shape,delay_params.scale,[sum(spikes) 1]));
            elseif i_outer==3
                delay = round(normrnd(delay_params.mean,delay_params.std,[sum(spikes) 1]));
            end
            spikes_est=find(responses(i_trial,:))-round(delay)';
            spikes_est=spikes_est( spikes_est>0);
            responses_MC(i_trial,spikes_est)=1;
            %sum(responses(i_trial,:))
        end
        
        % use the first event:
        if first_only
            responses_reg=responses_MC;
            responese_reg(:,:)=0;
            for i_trial = 1:size(responses_MC,1)
                spk_time= find(responses_MC(i_trial,:));
                if length(spk_time)>0
                    responese_reg(i_trial,spk_time(1))=1;
                end
            end
        else
            responses_reg = responses;
        end
        
        % Fit the lif-glm, and then
        in_params.g =   l23_cells_for_sim(i_template).g;
        % LIF-GLM fits
        %-------------------------------------%
        [stats_conv] = fit_lifglm_v3(responses_reg, stims,in_params,v_reset_known,first_only);
        
        % Output:
        MC_gain(i_MC)=stats_conv.beta(1);
    end
    cell_params.gain=median(MC_gain);
    
    [A,B,M_intensity_func]=Intensity_v7(outputM,... % whether to output M intensity
        stimuli_size, mpp,I_stimuli,... % data from exp
        cell_params, delay_params,... % estimated parameters
        funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
        n_stimuli_grid,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
        V_threshold,stimulus_threshold,first_only);
    estimated_intensities{i_outer}=M_intensity_func{1,1};
    
end


%%
%
figure(1)
colors=['r' 'g' 'b'];
t_grid = 1:length(I_e_vect);
for i_outer = 1:3
    plot(t_grid,estimated_intensities{i_outer},'col',colors(i_outer),'LineStyle','--','LineWidth',1)
    hold on;
    plot(t_grid,avg_spikes{i_outer},'col',colors(i_outer),'LineWidth',1)
    hold on;
end
xlabel('Time (1/20 ms)')
ylabel('Intensities')
hold off;
% saveas(1,'./Delay_check.png');

%% Iterative update of the LIF-GLM and the delays:
i_template=1;

for i_outer = 1:3
    responses = zeros(num_sample, length(I_e_vect));
    responses_true = zeros(num_sample, length(I_e_vect));
    
    stims = zeros(num_sample, length(I_e_vect));
    params_sim.gain =  l23_cells_for_sim(i_template).optical_gain;

    cell_params=params_sim;
    cell_params.gain_sd=0.01;
    
    if i_outer==1
        delay_params.type=1;
        delay_params.mean=0;
        delay_params.std=0.1;
    elseif i_outer==2
        delay_params.type=0;
        delay_params.shape=30*30/25;
        delay_params.scale=25/30;
    else
        delay_params.type=1;
        delay_params.mean=35;
        delay_params.std=15;
    end
    for i_trial = 1:num_sample
        k=stimuli_seq(i_trial);
        stim = I_e_vect*k;
        stims(i_trial,:) = stim;
        [V_vect, spikes]  = lif_glm_sim_v2(stim,params_sim,funcs);
        
        % shift the spikes give the delay
        
        if i_outer==1
            delay=0;
        elseif i_outer==2
            delay = round(gamrnd(delay_params.shape,delay_params.scale,[sum(spikes) 1]));
        elseif i_outer==3
            delay = round(normrnd(delay_params.mean,delay_params.std,[sum(spikes) 1]));
        end
        
        
        spikes_delayed=find(spikes)+round(delay)';
        spikes_delayed=spikes_delayed( spikes_delayed<time_max);
        responses(i_trial,spikes_delayed)=1;
        responses_true(i_trial,:)=spikes;
        
        %sum(responses(i_trial,:))
    end
    
    mpp=struct();
    for i=1:num_sample
        mpp(i).times=find(responses(i,:));
    end
    avg_spikes{i_outer}=mean(responses,1);
    
    cell_params.gain= mean([l23_cells_for_sim.optical_gain]);
    
    normalized_change = 1;
    convergence_epsilon=1e-4;
    maxit=100;
    num_iter=1;
    while (normalized_change > convergence_epsilon) & (num_iter < maxit)
        num_iter = num_iter+1;
        
        
        % Do not account for the delay when estimating the raw intensities 
      delay_params_est.type=1;
        delay_params_est.mean=0;
        delay_params_est.std=0.1;
        [A,B,M_intensity_func]=Intensity_v7(outputM,... % whether to output M intensity
            stimuli_size, mpp,I_stimuli,... % data from exp
            cell_params, delay_params_est,... % estimated parameters
            funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
            n_stimuli_grid,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
            V_threshold,stimulus_threshold,first_only);
        estimated_intensities_temp=M_intensity_func{1,1};
        
        responses_reg=zeros(num_sample, length(I_e_vect));
        for i_trial = 1:num_sample
            % first spike:
            spikes=find(responses(i_trial,:));
            if length(spikes)>0
                spike_first = spikes(1);
            end
            delay_seq =spike_first-(1:length(I_e_vect));
            if i_outer==1
                delay_prob = normpdf( delay_seq,delay_params.mean,delay_params.std);
            elseif i_outer==2
                delay_prob =gampdf(delay_seq,delay_params.shape,delay_params.scale);
            elseif i_outer==3
                delay_prob = normpdf(delay_seq,delay_params.mean,delay_params.std);
            end
            prod_prob = delay_prob.*estimated_intensities_temp';
            [~, spike_max]= max(prod_prob);
            responses_reg(i_trial,spike_max)=1;
            %sum(responses(i_trial,:))
        end
        % Update the delay distribution:
        delays=[];
         for i_trial = 1:num_sample
            % first spike:
            spikes=find(responses(i_trial,:));
            if length(spikes)>0
                spikes_adj = find(responses_reg(i_trial,:));
                delays =[delays spikes(1) - spikes_adj];
            end
         end
          
     if  delay_params.type == 1
       delay_params.mean=mean(delays);
        delay_params.std=std(delays);
    elseif delay_params.type == 0
        mean_temp = mean(delays);
        std_temp=std(delays);
        
        delay_params.scale = std_temp^2/mean_temp;
        delay_params.shape = mean_temp/delay_params.scale;
    end
        
        %
        % Fit the lif-glm, and then
        in_params.g =   l23_cells_for_sim(i_template).g;
        % LIF-GLM fits
        %-------------------------------------%
        [stats_conv] = fit_lifglm_v3(responses_reg, stims,in_params,v_reset_known,first_only);
        % Output:
        
        normalized_change = ((stats_conv.beta(1)-cell_params.gain)^2)/((cell_params.gain)^2);
        cell_params.gain=stats_conv.beta(1);
        fprintf('gain %d with %d changes\n',stats_conv.beta(1),normalized_change);
    end
      [A,B,M_intensity_func]=Intensity_v7(outputM,... % whether to output M intensity
            stimuli_size, mpp,I_stimuli,... % data from exp
            cell_params, delay_params,... % estimated parameters
            funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
            n_stimuli_grid,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
            stimulus_threshold,stimulus_threshold,first_only);
      
    estimated_intensities{i_outer}=M_intensity_func{1,1};
    avg_mle_spikes{i_outer}=mean(responses_reg,1);
    
end


%%
figure(1)
colors=['r' 'g' 'b'];
t_grid = 1:length(I_e_vect);
for i_outer = 1:3
    plot(t_grid,estimated_intensities{i_outer},'col',colors(i_outer),'LineStyle','--','LineWidth',1)
    hold on;
% plot(t_grid,estimated_intensities_temp,'col',colors(i_outer),'LineStyle','--','LineWidth',1)
%     hold on;

%     plot(t_grid,avg_mle_spikes{i_outer},'col',colors(i_outer),'LineWidth',1,'LineStyle',':')
%     hold on;
    plot(t_grid,avg_spikes{i_outer},'col',colors(i_outer),'LineWidth',1)
    hold on;
    
%     plot(t_grid,mean(responses_true,1),'col',colors(i_outer),'LineWidth',1,'LineStyle','-.')
%     hold on;
end
xlabel('Time (1/20 ms)')
ylabel('Intensities')
hold off;
 saveas(1,'./Delay_check.png');
