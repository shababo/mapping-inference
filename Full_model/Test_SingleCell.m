%% Testing the effect of delay on lif-glm & firing rate
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%%
v_th_known=15;
v_reset_known=-4000;
gain_initial=0.03;
num_per_power=100;
% delay_params.mean=35;
% delay_params.std=15;
delay_params.mean=0;
delay_params.std=0.1;

background_rate =4e-4;
% background_rate =4e-8;
gamma=1;
%%
rng(12242,'twister');
load('./Environments/l23_cells_for_sim.mat');
num_types_cell = length(l23_cells_for_sim);
% normalized the cell shapes
for i = 1:num_types_cell
    temp=l23_cells_for_sim(i).shape;
    temp_max = max(max(max(temp)));
    l23_cells_for_sim(i).shape = temp/temp_max;
end
i_template = 1;

params_sim.V_th= v_th_known;
params_sim.V_reset = v_reset_known;
params_sim.g =  l23_cells_for_sim(i_template).g;
params_sim.gain =  gain_initial;
funcs.invlink = @invlink_sig;%@(resp) log(1 + exp(resp));%@(x) exp(x);

%stimuli_seq = num_sample*(rand([num_sample 1])+0.5);
stimuli_size = [50*ones(num_per_power,1); 75*ones(num_per_power,1); 100*ones(num_per_power,1)];

%
load('./Environments/chrome-template-3ms.mat');
downsamp=1;
time_max= 300;
current_template=template(1:downsamp:time_max);

n_grid_voltage =1000;
voltage_threshold=-50;
dt=1;
n_grid_time = time_max;
t_factor=1;

sd_range=0.1;
n_delay_grid = 200;

stimulus_threshold=1e-4;
n_stimuli_grid=10;
gap_stimuli=0.05;

first_only=true;
V_threshold = -50;
n_trial = length(stimuli_size);

%%
responses = zeros(n_trial, length(current_template));
responses_raw = zeros(n_trial, length(current_template));
delay_raw = zeros(n_trial, length(current_template));
true_assignments = zeros(n_trial, length(current_template));
   
stims = zeros(n_trial, length(current_template));
cell_params=params_sim;
cell_params.gain_sd=0.01;

    mu_bg = 1/background_rate;
    spike_time_power=cell(4,1);
    spike_time=cell(3,1);
for i_trial = 1:n_trial
    k=stimuli_size(i_trial);
    stim = current_template*k;
    stims(i_trial,:) = stim;
    [V_vect, spikes]  =lif_glm_sim_v2(stim,params_sim,funcs);
    delay_vec=round(normrnd(delay_params.mean,delay_params.std,[sum(spikes) 1]));
    spikes_delay =find(spikes)+delay_vec;
    if spikes_delay < time_max
    if rand(1) < gamma
        responses(i_trial,spikes_delay)=1;
    true_assignments(i_trial,spikes_delay)=1;
    
    delay_raw(i_trial,find(spikes))=delay_vec;
    end
    end
    responses_raw(i_trial,:)=spikes;
%     if rand(1) < gamma
%     spike_time_power{floor((i_trial-1)/num_per_power)+1}=[spike_time_power{floor((i_trial-1)/num_per_power)+1} find( responses(i_trial,:))];
%     end
%     spike_time_power{floor((i_trial-1)/num_per_power)+1}=[spike_time_power{floor((i_trial-1)/num_per_power)+1} find(spikes)];
% %     
   
    % add background event: 
    R = exprnd(mu_bg);
    while R < time_max
        responses(i_trial, max(1,round(R)))=responses(i_trial, max(1,round(R)))+1;
        spike_time_power{4} =[spike_time_power{4} round(R)];
        R = R+exprnd(mu_bg);
        
    end
     
end
%% Find the minimal phi so that the cell fires at the highest power (direct stim)
 cell_params.gain = max(stimuli_size);
 gain_sequence = [1:100]/1e4; %from 0 to 0.01 
 n_stimuli_grid_temp=length(gain_sequence);
 
     [Stimuli_grid, Intensity_grid]=Intensity_v8(...
        gain_sequence,current_template,... % data from exp
        cell_params,... % estimated parameters
        funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
        n_stimuli_grid_temp,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
        V_threshold,stimulus_threshold,first_only);
 
    for i_grid = 1:n_stimuli_grid_temp
        if sum(Intensity_grid{i_grid}) > 0.99
            gain_lower_bound= gain_sequence(i_grid);
            break;
        end
    end
    

    %%
    normalized_change = 1;
    convergence_epsilon=1e-3;
    loss_type=2;
    initial_values = [gain_lower_bound  0.005 0.008 0.01 0.015 0.02 0.03 0.04 0.05];
    initial_values = [initial_values' 0.5*ones(length(initial_values),1)];
    % set a time threshold and throw away any spikes before this time point
    
    %      responses_MC=responses;
    %     stims_reg=stims;
    
    responses_reg=responses;
    stims_reg=stims;
    %
    in_params.g =   l23_cells_for_sim(i_template).g;
    % LIF-GLM fits
    %-------------------------------------%
    [stats_conv] = fit_lifglm_v5(responses_reg, stims_reg,in_params,...
        background_rate,v_reset_known,v_th_known,first_only,loss_type,...
        gain_lower_bound,initial_values);
    
    % Output:
    normalized_change = sqrt(((stats_conv(1)-cell_params.gain)^2)/((cell_params.gain)^2));
    cell_params.gain=stats_conv(1);
    
    [Stimuli_grid, Intensity_grid]=Intensity_v8(...
        stimuli_size,current_template,... % data from exp
        cell_params,... % estimated parameters
        funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
        n_stimuli_grid,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
        V_threshold,stimulus_threshold,first_only);

%%
stats_conv

    %%
xrange=gain_lower_bound + (1:400)/1e4;

  n_trial = size(stims,1);
    n_grid =size(stims,2);
    v_trace=zeros(n_trial,n_grid);
    for i_trial = 1:n_trial 
        v_trace(i_trial, 1) = 0;
        for i_t = 2:n_grid
            v_trace(i_trial, i_t) = sum( stims(i_trial, 1:(i_t-1)).*exp( ((1:(i_t-1)) -i_t+1)*in_params.g) );
        end
    end
    
    lklh=zeros(length(xrange),1);
        for i = 1:length(xrange)
  lklh(i)=lif_glm_firstspike_loglikelihood([xrange(i) 1], responses, v_trace,background_rate,v_th_known,linkfunc,loss_type);
        end
        %%
    figure(1)
    plot(xrange,lklh(:,1))
    xlim([min(xrange) max(xrange)])
    xlabel('Gain')
    ylabel('Log likelihood')
    title('True gain = 0.02')
    %xlim([0 0.003])
%     figure(2)
%     plot(xrange,log(lklh(:,2)))
%% Assignments:

soft_true=[];
true_time=[];

soft_bg=[];
bg_time =[];
for i_trial =1:n_trial
    events=find(responses(i_trial,:)>0);
    if length(events)>0
       for i = 1:length(events)
           if true_assignments(i_trial,events(i))==1
               soft_true=[soft_true soft_assignments(i_trial,events(i))];
               true_time=[true_time events(i)];   
           else
               soft_bg=[soft_bg soft_assignments(i_trial,events(i))];
               bg_time=[bg_time events(i)];
           end
       end
    end
end

plot(true_time,soft_true,'.','MarkerSize',10)
hold on;
plot(bg_time,soft_bg,'.','MarkerSize',10,'col','r')
hold off;
xlim([0,300])
ylim([0,1])
%%

soft_true=[];
true_time=[];

soft_bg=[];
bg_time =[];
for i_trial =1:n_trial
    events=find(responses_reg(i_trial,:)>0);
    if length(events)>0
       for i = 1:length(events)
           if true_assignments(i_trial,events(i))==1
               soft_true=[soft_true responses_reg(i_trial,events(i))];
               true_time=[true_time events(i)];   
           else
               soft_bg=[soft_bg responses_reg(i_trial,events(i))];
               bg_time=[bg_time events(i)];
           end
       end
    end
end

plot(true_time,soft_true,'.','MarkerSize',10)
hold on;
plot(bg_time,soft_bg,'.','MarkerSize',10,'col','r')
hold off;
xlim([0,300])
ylim([0,1])

%% Draw the spike times 
figure(2)
plot(sort(spike_time_power{1}),1:length(spike_time_power{1}),...
    '.','col','g','MarkerSize',20)
hold on;
plot(sort(spike_time_power{2}),1:length(spike_time_power{2}),...
    '.','col','b','MarkerSize',20)
hold on;
plot(sort(spike_time_power{3}),1:length(spike_time_power{3}),...
'.','col','r','MarkerSize',20)
hold off;
ylim([0,num_per_power]);
xlim([0,300]);
line([mean(spike_time_power{1}) mean(spike_time_power{1}) ], [0 20],'col','g')
line([mean(spike_time_power{2}) mean(spike_time_power{2}) ], [0 20],'col','b')
line([mean(spike_time_power{3}) mean(spike_time_power{3}) ], [0 20],'col','r')
xlabel('Time (1/20 ms)');
ylabel('Ranks (irrelevant)');

% saveas(2,'./sorted_time.png');
% saveas(2,'./sorted_time_sim.png');

%%
figure(1)
colors=['g' 'b' 'r'];
t_grid = 1:length(current_template);
for i_outer = 1:3
    avg_spikes{i_outer}=mean(responses( (i_outer-1)*num_per_power +(1:num_per_power),: ));
    %     plot(t_grid,Intensity_grid{i_outer},'col',colors(i_outer),'LineStyle','--','LineWidth',1)
    %     hold on;
    plot(t_grid,avg_spikes{i_outer},'col',colors(i_outer),'LineStyle','--','LineWidth',1)
    hold on;
     avg_spikes_true{i_outer}=mean(responses_raw( (i_outer-1)*num_per_power +(1:num_per_power),: ));
    plot(t_grid,avg_spikes_true{i_outer},'col',colors(i_outer),'LineStyle',':','LineWidth',1)
    hold on;
   
end
 xlim([0,300]);
    xlabel('Time (1/20 ms)')
    ylabel('Intensities')
    hold off;
%%

figure(1)
colors=['g' 'b' 'r'];
t_grid = 1:length(current_template);
for i_outer = 1:3
    figure(i_outer)
%     avg_spikes{i_outer}=mean(responses_reg( (i_outer-1)*num_per_power +(1:num_per_power),: ));
    avg_spikes_true{i_outer}=mean(responses_raw( (i_outer-1)*num_per_power +(1:num_per_power),: ));
    
    plot(t_grid,Intensity_grid{i_outer},'col',colors(i_outer),'LineWidth',1)
    hold on;
%     plot(t_grid,avg_spikes{i_outer},'col',colors(i_outer),'LineStyle','--','LineWidth',1)
%     hold on;
    plot(t_grid,avg_spikes_true{i_outer},'col',colors(i_outer),'LineStyle',':','LineWidth',1)
    hold on;
    
    xlim([0,300]);
    xlabel('Time (1/20 ms)')
    ylabel('Intensities')
    hold off;
end
% saveas(3,'./Delay_bg.png');

%%

figure(1)
colors=['r' 'g' 'b'];
t_grid = 1:length(current_template);
for i_outer = 1:3
      plot(t_grid,Intensity_grid{i_outer},'col',colors(i_outer),'LineStyle','--','LineWidth',1)
    hold on;
  
    avg_spikes_true{i_outer}=mean(responses_raw( (i_outer-1)*num_per_power +(1:num_per_power),: ));
    
    plot(t_grid,avg_spikes_true{i_outer},'col',colors(i_outer),'LineStyle',':','LineWidth',1)
    hold on;
   
end
xlim([0,300]);
xlabel('Time (1/20 ms)')
    
%% Iterative update of the LIF-GLM and the delays:
% i_template=1;
% 
% for i_outer = 1:3
%     responses = zeros(num_sample, length(I_e_vect));
%     responses_true = zeros(num_sample, length(I_e_vect));
%     
%     stims = zeros(num_sample, length(I_e_vect));
%     params_sim.gain =  l23_cells_for_sim(i_template).optical_gain;
%     
%     cell_params=params_sim;
%     cell_params.gain_sd=0.01;
%     
%     if i_outer==1
%         delay_params.type=1;
%         delay_params.mean=0;
%         delay_params.std=0.1;
%     elseif i_outer==2
%         delay_params.type=0;
%         delay_params.shape=30*30/25;
%         delay_params.scale=25/30;
%     else
%         delay_params.type=1;
%         delay_params.mean=35;
%         delay_params.std=15;
%     end
%     for i_trial = 1:num_sample
%         k=stimuli_seq(i_trial);
%         stim = I_e_vect*k;
%         stims(i_trial,:) = stim;
%         [V_vect, spikes]  = lif_glm_sim_v2(stim,params_sim,funcs);
%         
%         % shift the spikes give the delay
%         
%         if i_outer==1
%             delay=0;
%         elseif i_outer==2
%             delay = round(gamrnd(delay_params.shape,delay_params.scale,[sum(spikes) 1]));
%         elseif i_outer==3
%             delay = round(normrnd(delay_params.mean,delay_params.std,[sum(spikes) 1]));
%         end
%         
%         
%         spikes_delayed=find(spikes)+round(delay)';
%         spikes_delayed=spikes_delayed( spikes_delayed<time_max);
%         responses(i_trial,spikes_delayed)=1;
%         responses_true(i_trial,:)=spikes;
%         
%         %sum(responses(i_trial,:))
%     end
%     
%     mpp=struct();
%     for i=1:num_sample
%         mpp(i).times=find(responses(i,:));
%     end
%     avg_spikes{i_outer}=mean(responses,1);
%     
%     cell_params.gain= mean([l23_cells_for_sim.optical_gain]);
%     
%     normalized_change = 1;
%     convergence_epsilon=1e-4;
%     maxit=100;
%     num_iter=1;
%     while (normalized_change > convergence_epsilon) & (num_iter < maxit)
%         num_iter = num_iter+1;
%         
%         
%         % Do not account for the delay when estimating the raw intensities
%         delay_params_est.type=1;
%         delay_params_est.mean=0;
%         delay_params_est.std=0.1;
%         [A,B,M_intensity_func]=Intensity_v7(outputM,... % whether to output M intensity
%             stimuli_size, mpp,I_stimuli,... % data from exp
%             cell_params, delay_params_est,... % estimated parameters
%             funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
%             n_stimuli_grid,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
%             V_threshold,stimulus_threshold,first_only);
%         estimated_intensities_temp=M_intensity_func{1,1};
%         
%         responses_reg=zeros(num_sample, length(I_e_vect));
%         for i_trial = 1:num_sample
%             % first spike:
%             spikes=find(responses(i_trial,:));
%             if length(spikes)>0
%                 spike_first = spikes(1);
%             end
%             delay_seq =spike_first-(1:length(I_e_vect));
%             if i_outer==1
%                 delay_prob = normpdf( delay_seq,delay_params.mean,delay_params.std);
%             elseif i_outer==2
%                 delay_prob =gampdf(delay_seq,delay_params.shape,delay_params.scale);
%             elseif i_outer==3
%                 delay_prob = normpdf(delay_seq,delay_params.mean,delay_params.std);
%             end
%             prod_prob = delay_prob.*estimated_intensities_temp';
%             [~, spike_max]= max(prod_prob);
%             responses_reg(i_trial,spike_max)=1;
%             %sum(responses(i_trial,:))
%         end
%         % Update the delay distribution:
%         delays=[];
%         for i_trial = 1:num_sample
%             % first spike:
%             spikes=find(responses(i_trial,:));
%             if length(spikes)>0
%                 spikes_adj = find(responses_reg(i_trial,:));
%                 delays =[delays spikes(1) - spikes_adj];
%             end
%         end
%         
%         if  delay_params.type == 1
%             delay_params.mean=mean(delays);
%             delay_params.std=std(delays);
%         elseif delay_params.type == 0
%             mean_temp = mean(delays);
%             std_temp=std(delays);
%             
%             delay_params.scale = std_temp^2/mean_temp;
%             delay_params.shape = mean_temp/delay_params.scale;
%         end
%         
%         %
%         % Fit the lif-glm, and then
%         in_params.g =   l23_cells_for_sim(i_template).g;
%         % LIF-GLM fits
%         %-------------------------------------%
%         [stats_conv] = fit_lifglm_v3(responses_reg, stims,in_params,v_reset_known,first_only);
%         % Output:
%         
%         normalized_change = ((stats_conv.beta(1)-cell_params.gain)^2)/((cell_params.gain)^2);
%         cell_params.gain=stats_conv.beta(1);
%         fprintf('gain %d with %d changes\n',stats_conv.beta(1),normalized_change);
%     end
%     [A,B,M_intensity_func]=Intensity_v7(outputM,... % whether to output M intensity
%         stimuli_size, mpp,I_stimuli,... % data from exp
%         cell_params, delay_params,... % estimated parameters
%         funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
%         n_stimuli_grid,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
%         stimulus_threshold,stimulus_threshold,first_only);
%     
%     estimated_intensities{i_outer}=M_intensity_func{1,1};
%     avg_mle_spikes{i_outer}=mean(responses_reg,1);
%     
% end
% 
% 
% %%
% figure(1)
% colors=['r' 'g' 'b'];
% t_grid = 1:length(I_e_vect);
% for i_outer = 1:3
%     plot(t_grid,estimated_intensities{i_outer},'col',colors(i_outer),'LineStyle','--','LineWidth',1)
%     hold on;
%     % plot(t_grid,estimated_intensities_temp,'col',colors(i_outer),'LineStyle','--','LineWidth',1)
%     %     hold on;
%     
%     %     plot(t_grid,avg_mle_spikes{i_outer},'col',colors(i_outer),'LineWidth',1,'LineStyle',':')
%     %     hold on;
%     plot(t_grid,avg_spikes{i_outer},'col',colors(i_outer),'LineWidth',1)
%     hold on;
%     
%     %     plot(t_grid,mean(responses_true,1),'col',colors(i_outer),'LineWidth',1,'LineStyle','-.')
%     %     hold on;
% end
% xlabel('Time (1/20 ms)')
% ylabel('Intensities')
% hold off;
% saveas(1,'./Delay_check.png');
% % 

%% Waste code on soft assignments 
   
  %     nonbackground_prob=zeros(n_trial, length(current_template));
%     for i_trial = 1:n_trial
% %         if num_iter == 2
% %         estimated_intensities_temp=ones(length(current_template),1);
% %         else
%             [~,min_idx]=min(abs(Stimuli_grid-stimuli_size(i_trial)*cell_params.gain));
%         estimated_intensities_temp=Intensity_grid{min_idx};
% %         end
%         spikes=find(responses(i_trial,:));
%         if length(spikes)>0
%             for i = 1:length(spikes)
%                 spike_time = spikes(i);
%                 delay_seq =spike_time-(1:length(current_template));
%                 delay_prob = normpdf(delay_seq,delay_params.mean,delay_params.std);
%                 prod_prob = delay_prob.*estimated_intensities_temp';
%                 nonbackground_prob(i_trial, spike_time)= sum(prod_prob);
%             end
%         end
%     end
  % Calculate the probability of spikes
% %     
% %     % Calculate the soft assignments
%     soft_assignments=zeros(n_trial, length(current_template));
%     for i_trial = 1:n_trial
%         spikes=find(responses(i_trial,:));
%         if length(spikes)>0
%             for i = 1:length(spikes)
%                 spike_time = spikes(i);
%                 soft_assignments(i_trial, spike_time)= gamma*nonbackground_prob(i_trial, spike_time)/(nonbackground_prob(i_trial, spike_time)+background_rate);
%             end
%         end
%     end


%     
%     
    % Draw one MC sample based on the soft assignment 
%     responses_MC=zeros(n_trial, length(current_template));
%      for i_trial = 1:n_trial
%         potential_spikes=find(soft_assignments(i_trial,:)>0);
%         if length(potential_spikes)>0
%              for i = 1:length(potential_spikes)
%                 spike_indicator=rand(1) < soft_assignments(i_trial, potential_spikes(i));
%                 responses_MC(i_trial,potential_spikes(i))= spike_indicator;
%             end
%         end
%      end 


% responses_MC=soft_assignments;
% % responses_reg=responses_MC;
% stims_reg=stims;
%     responses_reg=zeros(n_trial, length(current_template));
%     delays=[];
%     counter=1;
%     for i_trial = 1:n_trial
%         
%         [~,min_idx]=min(abs(Stimuli_grid-stimuli_size(i_trial)*cell_params.gain));
%         estimated_intensities_temp=Intensity_grid{min_idx};
%         % first spike:
%           spikes=find(responses_MC(i_trial,:)>0);
%       
%         if length(spikes)>0
%             for i= 1:length(spikes)
%             spike_first = spikes(i);
%             delay_seq =spike_first-(1:length(current_template));
%             delay_prob = normpdf(delay_seq,delay_params.mean,delay_params.std);
%             prod_prob = delay_prob.*estimated_intensities_temp';
%             [~, spike_max]= max(prod_prob);
%             delays =[delays spikes(i) - spike_max];
%             responses_reg(i_trial,spike_max)=responses_MC(i_trial,spike_first);
%             %counter=counter+1;
%             end
%         end
%            
%     end
    
    
