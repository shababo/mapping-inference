addpath(genpath('../../../mapping-inference'));
%% Parameters:
d=10; % distance between the two cells
%% Cell locations
cell_locations=zeros(2,3);
cell_locations(2,1)=d;
n_cell=2;
%% Cell parameters
num_sim=4;
num_seed=20;
num_threshold = 2;
for i_setting = 1:num_sim
    background_rate=1e-4;
    v_th_known=15*ones([n_cell,1]);v_reset_known=-1e4*ones([n_cell,1]);
    g_truth = 0.02*ones([n_cell,1]);
    funcs.invlink = @invlink_sig;%@(resp) log(1 + exp(resp));%@(x) exp(x);
    
        %     gain_truth=  gain_bound.low+ (gain_bound.up-gain_bound.low)*rand([1 2]);
        gain_truth=[0.015; 0.02];
        % Draw gains
        gamma_truth=zeros(1,2);
        gamma_truth(1)= i_setting/10;
    
        num_trace_back=3;
    for i_seed = 1:num_seed
        rng(i_seed,'twister');
        %%  Load the current template
        load('../Environments/chrome-template-3ms.mat');
        downsamp=1;max_time=300;power_level = 30:10:100;
        num_power_level=length(power_level);
        current_template=template(1:downsamp:max_time);
        t_vect= 1:1:max_time;
        time_max=max_time;
        %% Preprocessing
        % Calculate the firing probability
        delay_params.type=2; %1: normal; 2: gamma
        delay_params.mean=58; delay_params.std=15;
        delay_params.delayed=true; delay_params.n_grid=200;
        cell_params.gain_template = 1; % for calculation
        
        cell_params.g=0.02;cell_params.v_th_known=15;
        linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};
        gain_template=0.02;
        stim_scale=4/gain_template;
        stim_grid = (1:1000)/stim_scale;
        stim_unique=(1:1000)/stim_scale/gain_template;
        [prob_trace_full,v_trace_full] = get_first_spike_intensity(...
            linkfunc,...
            current_template,stim_grid,cell_params,delay_params);
        prob_trace=sum(prob_trace_full,2);
        eff_stim_threshold=stim_grid(min(find(prob_trace>0.01)));
        fire_stim_threshold=stim_grid(min(find(prob_trace>0.99)));
        
        %% Select the stimulation locations
        % Load the shape template
        load('../Environments/l23_template_cell.mat');
        %l23_average_shape
        temp=l23_average_shape;temp_max = max(max(max(temp)));
        l23_average_shape = temp/temp_max;shape_template=l23_average_shape;
        r1=5;r2=10;r3=15;num_per_grid=12;
        num_per_grid_dense=16;
        % The list of all related cells :
        grid_jitters = zeros(num_per_grid,2);
        for i_grid = 1:num_per_grid
            grid_jitters(i_grid,:)=[sin(2*pi*i_grid/num_per_grid) cos(2*pi*i_grid/num_per_grid)];
        end
        grid_jitters=[grid_jitters zeros(num_per_grid,1)];
        
        % Calculate the stimulation locations
        target_locations = zeros(n_cell*(2*num_per_grid+1),3);
        for i_cell=1:n_cell
            nucleus_loc=cell_locations(i_cell,:);
            grid_locs=nucleus_loc;
            grid_locs=[grid_locs;bsxfun(@plus,grid_jitters*r1,nucleus_loc)];
            grid_locs=[grid_locs;bsxfun(@plus,grid_jitters*r2,nucleus_loc)];
            target_idx=(i_cell-1)*(2*num_per_grid+1) +(1: (2*num_per_grid+1));
            target_locations(target_idx,:) = grid_locs;
        end
        % target_locations(:,3)= mean(cell_locations(target_cell_list.primary,3));
        
        %plot(target_locations{this_plane}(:,2),target_locations{this_plane}(:,1),'.')
        
        
        cell_params.locations =  cell_locations;
        cell_params.shape_gain = ones(n_cell,1);
        cell_template = struct();
        cell_template.shape= shape_template;
        % [pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
        %     cell_template,target_locations);
        [pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
            cell_template,target_locations);
        %power_selected = power_level(1)*ones([size(target_locations,1) 1]);
        loc_to_cell = zeros(size(pi_target,2),1);
        loc_to_cell(1:(2*num_per_grid+1))=1;
        loc_to_cell((2*num_per_grid+1)+ (1:(2*num_per_grid+1)))=2;
        power_selected=zeros(n_cell*(2*num_per_grid+1),1);
        power_sd=zeros(n_cell*(2*num_per_grid+1),1);
        
        %% Parameters in the design stage
        % Design parameters
        n_spots_per_trial = 1;
        % Need to run sims to check how these parameters affect the results
        % Prior distribution
        prior_pi0=0.8;
        
        % Initialize tuning parameters in the VI
        maxit=1000;
        S=200;epsilon=0.01;eta_logit=1;eta_beta=1;eta_max=2;
        background_rt=background_rate*time_max; % raw probability of firing within a trial
        
        gamma_estimates = 0.5*ones(n_cell,1);% for drawing samples (not really used)
        prob_weight=0;
        id_continue=1;% an indicator
        
        % lklh_func=@calculate_likelihood_sum_bernoulli; % likelihood function is
        % specificed when fitting the working model
        
        stim_threshold = 10;
        % bounds of the gamma:
        gain_bound.up=0.03;gain_bound.low=0.005;
        
        
        % Initialize storage
        mean_gamma_current=zeros(n_cell,1);
        mean_gain_current=gain_template*ones(n_cell,1);
        gamma_path=zeros(n_cell,1);
        gain_path=zeros(n_cell,1);
        var_gamma_path=zeros(n_cell,1);
        
        %% Simulate the prior distribution of gain:
        
        var_gain_prior=1;gain_bias=0.005;K=5;
        n_replicates=1;cell_list= 1:2; trial_max=100;
        iter=1;n_MC_samples=100;epislon=10;
        connected_threshold=0.4;
%         %%
%         % bounds of the gamma:
%         gain_bound.up=0.03;gain_bound.low=0.005;
%         var_alpha_initial=1;var_beta_initial=1.78;
%         var_alpha_gain_initial=log( (gain_truth+gain_bias*(2*(rand(1)-0.5))- gain_bound.low)./(gain_bound.up-gain_truth));
%         var_beta_gain_initial=var_gain_prior; % uncertainty of the prior
%         % The prioir info serves as the first variational distributions
%         variational_params_path.alpha=var_alpha_initial*ones(n_cell,1);
%         variational_params_path.beta=var_beta_initial*ones(n_cell,1);
%         variational_params_path.alpha_gain=var_alpha_gain_initial;
%         variational_params_path.beta_gain=var_beta_gain_initial*ones(n_cell,1);
%         run('Twoneurons_naive.m');
%         n_events_naive
        %%
        K=5;
    C_threshold = 0.01;
        % bounds of the gamma:
        gain_bound.up=0.03;gain_bound.low=0.005;
        var_alpha_initial=1;var_beta_initial=1.78;
        var_alpha_gain_initial=log( (gain_truth+gain_bias*(2*(rand(1)-0.5)) - gain_bound.low)./(gain_bound.up-gain_truth));
        var_beta_gain_initial=var_gain_prior; % uncertainty of the prior
        % The prioir info serves as the first variational distributions
        variational_params_path.alpha=var_alpha_initial*ones(n_cell,1);
        variational_params_path.beta=var_beta_initial*ones(n_cell,1);
        variational_params_path.alpha_gain=var_alpha_gain_initial;
        variational_params_path.beta_gain=var_beta_gain_initial*ones(n_cell,1);
        run('Twoneurons_optimal.m');
        n_events_optimal
%         %% Use the full model:
%         
%         % bounds of the gamma:
%         gain_bound.up=0.03;gain_bound.low=0.005;
%         var_alpha_initial=1;var_beta_initial=1.78;
%         var_alpha_gain_initial=log( (gain_truth+gain_bias*(2*(rand(1)-0.5)) - gain_bound.low)./(gain_bound.up-gain_truth));
%         var_beta_gain_initial=var_gain_prior; % uncertainty of the prior
%         % The prioir info serves as the first variational distributions
%         variational_params_path.alpha=var_alpha_initial*ones(n_cell,1);
%         variational_params_path.beta=var_beta_initial*ones(n_cell,1);
%         variational_params_path.alpha_gain=var_alpha_gain_initial;
%         variational_params_path.beta_gain=var_beta_gain_initial*ones(n_cell,1);
%         run('Twoneurons_fullmodel.m');
%         n_events_full
        
        %% Use spike and slab distribution 
        K=5;
params_path_spike=struct([]);
        for i_threshold= 1:num_threshold
               mean_gamma_current=zeros(n_cell,1);
        mean_gain_current=gain_template*ones(n_cell,1);
        gamma_path=zeros(n_cell,1);
        gain_path=zeros(n_cell,1);
        var_gamma_path=zeros(n_cell,1);
     
    C_threshold = 0.05*i_threshold;
        % bounds of the gamma:
        gain_bound.up=0.03;gain_bound.low=0.005;
        var_alpha_initial=0;var_beta_initial=1.78;
        var_alpha_gain_initial=log( (gain_truth+gain_bias*(2*(rand(1)-0.5)) - gain_bound.low)./(gain_bound.up-gain_truth));
        var_beta_gain_initial=var_gain_prior; % uncertainty of the prior
        % The prioir info serves as the first variational distributions
      
    variational_params_path.pi=0.5*ones(n_cell,1);
    variational_params_path.p_logit=zeros(n_cell,1);
     variational_params_path.alpha=var_alpha_initial*ones(n_cell,1);
        variational_params_path.beta=var_beta_initial*ones(n_cell,1);
        variational_params_path.alpha_gain=var_alpha_gain_initial;
        variational_params_path.beta_gain=var_beta_gain_initial*ones(n_cell,1);
        
          variational_params_path_full.pi=0.5*ones(n_cell,1);
    variational_params_path_full.p_logit=zeros(n_cell,1);
     variational_params_path_full.alpha=var_alpha_initial*ones(n_cell,1);
        variational_params_path_full.beta=var_beta_initial*ones(n_cell,1);
        variational_params_path_full.alpha_gain=var_alpha_gain_initial;
        variational_params_path_full.beta_gain=var_beta_gain_initial*ones(n_cell,1);
        run('Twoneurons_spike.m');
    
        

params_path_spike(i_threshold).pi=variational_params_path.pi;
params_path_spike(i_threshold).alpha=variational_params_path.alpha;
params_path_spike(i_threshold).beta = variational_params_path.beta;
params_path_spike(i_threshold).alpha_gain= variational_params_path.alpha_gain;
params_path_spike(i_threshold).beta_gain= variational_params_path.beta_gain;

pi_path_spike{i_threshold}=variational_params_path.pi;
var_gamma_path_spike{i_threshold}=var_gamma_path;
var_gain_path_spike{i_threshold}=var_gain_path;
gamma_path_spike{i_threshold}=gamma_path;
gain_path_spike{i_threshold}=gain_path;
n_events_spike{i_threshold}=n_events;
     

params_path_spike(num_threshold+i_threshold).pi=variational_params_path_full.pi;
params_path_spike(num_threshold+i_threshold).alpha=variational_params_path_full.alpha;
params_path_spike(num_threshold+i_threshold).beta = variational_params_path_full.beta;
params_path_spike(num_threshold+i_threshold).alpha_gain= variational_params_path_full.alpha_gain;
params_path_spike(num_threshold+i_threshold).beta_gain= variational_params_path_full.beta_gain;

pi_path_spike{num_threshold+i_threshold}=variational_params_path_full.pi;
var_gamma_path_spike{num_threshold+i_threshold}=var_gamma_path_full;
var_gain_path_spike{num_threshold+i_threshold}=var_gain_path_full;
gamma_path_spike{num_threshold+i_threshold}=gamma_path_full;
gain_path_spike{num_threshold+i_threshold}=gain_path_full;

n_events_spike{num_threshold+i_threshold}=n_events;
        end
        %%
        save(strcat('./matfiles/Oct18/','Sim', num2str(i_setting),'Seed',num2str(i_seed),'.mat'),...
            'gamma_path_optimal','var_gamma_path_optimal',...
            'gain_path_optimal','var_gain_path_optimal',...
            'params_path_optimal','fitting_time_optimal',...
            'mpp_optimal','trials_locations_optimal','trials_powers_optimal',...
            'gamma_path_spike','var_gamma_path_spike',...
            'gain_path_spike','var_gain_path_spike','pi_path_spike',...
            'params_path_spike','fitting_time_spike',...
            'mpp_spike','trials_locations_spike','trials_powers_spike',...
            'K','trial_max','n_cell','gain_truth','gamma_truth')
        
    end
end
%%
num_sim=4;
num_seed=20;
n_cell=2;
K=5;
trial_max=100;
num_threshold=2;

gamma_path_all=zeros(num_sim,num_seed,7,n_cell,(trial_max/(K *n_cell)));
gain_path_all=zeros(num_sim,num_seed,7,n_cell,(trial_max/(K *n_cell)));
var_gamma_path_all=zeros(num_sim,num_seed,7,n_cell,(trial_max/(K *n_cell)));
var_gain_path_all=zeros(num_sim,num_seed,7,n_cell,(trial_max/(K *n_cell)));
fitting_time_all=zeros(num_sim,num_seed,7,(trial_max/(K *n_cell)));
nonzero_prob=zeros(num_sim,num_seed,6,n_cell,(trial_max/(K *n_cell)));
for i_setting = 1:num_sim
      %     gain_truth=  gain_bound.low+ (gain_bound.up-gain_bound.low)*rand([1 2]);
        gain_truth=[0.015; 0.02];
        % Draw gains
        gamma_truth=zeros(1,2);
        gamma_truth(1)= i_setting/10;
      
    for i_seed = 1:num_seed
        % i_sim=0;i_seed=1;
        
        % Load data
        load(strcat('./matfiles/Oct18/','Sim', num2str(i_setting),'Seed',num2str(i_seed),'.mat'))
            
        gamma_path_all(i_setting,i_seed,2*num_threshold +1,:,:)=(gamma_path_optimal(:,2:end)-gamma_truth').^2;
        gain_path_all(i_setting,i_seed,2*num_threshold +1,:,:)=(gain_path_optimal(:,2:end)-gain_truth).^2;
        var_gamma_path_all(i_setting,i_seed,2*num_threshold +1,:,:)=var_gamma_path_optimal(:,2:end);
        var_gain_path_all(i_setting,i_seed,2*num_threshold +1,:,:)=var_gain_path_optimal(:,2:end);
      fitting_time_all(i_setting,i_seed,2*num_threshold +1,:)=fitting_time_optimal;
      
      for i_threshold=1:num_threshold 
        gamma_path_all(i_setting,i_seed,i_threshold,:,:)=(gamma_path_spike{i_threshold}(:,2:end)-gamma_truth').^2;
        gain_path_all(i_setting,i_seed,i_threshold,:,:)=(gain_path_spike{i_threshold}(:,2:end)-gain_truth).^2;
        var_gamma_path_all(i_setting,i_seed,i_threshold,:,:)=var_gamma_path_spike{i_threshold}(:,2:end);
        var_gain_path_all(i_setting,i_seed,i_threshold,:,:)=var_gain_path_spike{i_threshold}(:,2:end);
        nonzero_prob(i_setting,i_seed,i_threshold,:,:)=(1-pi_path_spike{i_threshold}(:,2:end));
      end

      num_batch=size(pi_path_spike{num_threshold+i_threshold},2)-1;
      for i_threshold=1:num_threshold 
        gamma_path_all(i_setting,i_seed,num_threshold+i_threshold,:,1:num_batch)=(gamma_path_spike{num_threshold+i_threshold}(:,2:11)-gamma_truth').^2;
        gain_path_all(i_setting,i_seed,num_threshold+i_threshold,:,1:num_batch)=(gain_path_spike{num_threshold+i_threshold}(:,2:end)-gain_truth).^2;
        var_gamma_path_all(i_setting,i_seed,num_threshold+i_threshold,:,1:num_batch)=var_gamma_path_spike{num_threshold+i_threshold}(:,2:end);
        var_gain_path_all(i_setting,i_seed,num_threshold+i_threshold,:,1:num_batch)=var_gain_path_spike{num_threshold+i_threshold}(:,2:end);
        nonzero_prob(i_setting,i_seed,num_threshold+i_threshold,:,1:num_batch)=(1-pi_path_spike{num_threshold+i_threshold}(:,2:end));
      end
   end
end
%% Take average

mean_gamma_path_all=mean(gamma_path_all,2);
% mean_var_gamma_all=mean(var_gamma_path_all,2);
mean_gain_path_all=mean(gain_path_all,2);
% mean_var_gain_all=mean(var_gain_path_all,2);
mean_fitting_time_all=mean(fitting_time_all,2);
mean_nonzero_prob=mean(nonzero_prob,2);
%% Plot:
K=5;
K2=5;
num_threshold=2;
line_styles={'-' ':'};
color_list=zeros(4,3,4);
color_list(1,:,:)=[ [1 0 0 0.3];[1 0 0 0.3];[1 0 0 0.3]];
color_list(2,:,:)=[ [0 1 0 0.3];[0 1 0 0.3];[0 1 0 0.3]];
color_list(3,:,:)=[ [1 0 0 0.7];[1 0 0 0.3];[1 0 0 0.3]];
color_list(4,:,:)=[ [0 1 0 0.7];[0 1 0 0.3];[0 1 0 0.3]];

% color_list(4,:,:)=[ [1 1 0 0.7];[1 1 0 0.3];[1 1 0 0.3]];

for i_setting = 1:num_sim
       gain_truth=[0.015; 0.02];
        % Draw gains
        gamma_truth=zeros(1,2);
        gamma_truth(1)= i_setting/10;
     
    
    % i_sim=0;i_seed=1;
%     
    figure(1)
    hold on;
    for i_cell = 1:n_cell
        for k=1:num_threshold
            line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                reshape(mean_gamma_path_all(i_setting,1,k,i_cell,:) ,[],1),...
                'Color',color_list(k,1,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
        end
        for k=1:num_threshold
            line((1: (trial_max/(K*n_cell)))*K*n_cell,...
                reshape(mean_gamma_path_all(i_setting,1,num_threshold+k,i_cell,:) ,[],1),...
                'Color',color_list(k+num_threshold,1,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
        end
        
    end
   
%     line([K*n_cell trial_max], [gamma_truth(1) gamma_truth(1)],...
%         'Color',[0 0 0 0.5],...
%         'LineStyle','-','LineWidth',3)
%     line([K*n_cell trial_max], [gamma_truth(2) gamma_truth(2)],...
%         'Color',[0 0 0 0.5],...
%         'LineStyle',':','LineWidth',3)
    
    ylim([0 0.2]);
      title( strcat('\gamma_1 =', num2str(gamma_truth(1)),', \gamma_2=0.'),'FontSize',15); 
   
    xlabel('Number of trials','FontSize',15);
    ylabel('E||E(\gamma|data)-\gamma||^2','FontSize',25);
    hold off;
%     %
%     
%     
%     figure(2)
%     hold on;
%     for i_cell = 1:n_cell
%         for k=1:num_models
%             line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%                 reshape(mean_gain_path_all(i_setting,1,k,i_cell,:) ,[],1),...
%                 'Color',color_list(k,1,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
%         end
%     end
% %     line([K*n_cell trial_max], [gain_truth(1) gain_truth(1)],...
% %         'Color',[0 0 0 0.5],...
% %         'LineStyle','-','LineWidth',3)
%     
% %     line([K*n_cell trial_max], [gain_truth(2) gain_truth(2)],...
% %         'Color',[0 1 0 0.5],...
% %         'LineStyle',':','LineWidth',3)
%     ylim([0 2e-5]);
%      title( strcat('\gamma_1 =', num2str(gamma_truth(1)),', \gamma_2=0.'),'FontSize',15); 
%    
%     xlabel('Number of trials','FontSize',15);
%     ylabel('E||E(\phi|data)-\phi||^2','FontSize',25);
%     hold off;
%     
%     figure(3)
%     hold on;
%     for i_cell = 1:n_cell
%         for k=1:num_models
%             line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%                 reshape(mean_var_gamma_all(i_setting,1,k,i_cell,:) ,[],1),...
%                 'Color',color_list(k,1,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
%         end
%     end
%     ylim([0 5e-2]);
%      title( strcat('\gamma_1 =', num2str(gamma_truth(1)),', \gamma_2=0.'),'FontSize',15); 
%    
%     xlabel('Number of trials','FontSize',15);
%     ylabel('var[\gamma|data]','FontSize',25);
%     hold off;
%     
%     
%     figure(4)
%     hold on;
%     for i_cell = 1:n_cell
%         for k=1:num_models
%             line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%                 reshape(mean_var_gain_all(i_setting,1,k,i_cell,:) ,[],1),...
%                 'Color',color_list(k,1,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
%         end
%     end
%     ylim([0 5e-5]);
%     
%     xlabel('Number of trials','FontSize',15);
%     ylabel('var[\phi|data]','FontSize',25);
%    title( strcat('\gamma_1 =', num2str(gamma_truth(1)),', \gamma_2=0.'),'FontSize',15); 
%      hold off;
    
     
%     figure(5)
%     hold on;
%     for k=1:num_models
%         line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%             reshape(mean_fitting_time_all(i_setting,1,k,:),[],1),...
%             'Color',color_list(k,1,:),'LineStyle','-','LineWidth',3)
%     end
%     ylim([0 4]);
%     
%     xlabel('Number of trials','FontSize',15);
%     ylabel('Comp. time','FontSize',25);
%     title( strcat('\gamma_1 =', num2str(gamma_truth(1)),', \gamma_2=0.'),'FontSize',15); 
%     hold off;
    
      
    figure(6)
    hold on;
    for i_cell = 1:n_cell
%            for k=1:num_threshold
%      
%             line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%                 reshape(mean_nonzero_prob(i_setting,1,k,i_cell,:) ,[],1),...
%                 'Color',color_list(k,1,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
%            end
%            for k=1:num_threshold
%             line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%                 reshape(mean_nonzero_prob(i_setting,1,num_threshold+k,i_cell,1:num_batch) ,[],1),...
%                 'Color',color_list(num_threshold+k,1,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
%            end
           for i_sim= 1:num_seed
           for k=1:num_threshold
            line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                reshape(nonzero_prob(i_setting,i_sim,k,i_cell,:) ,[],1),...
                'Color',color_list(k,1,:),'LineStyle',line_styles{i_cell},'LineWidth',1)
           end
%            for k=1:num_threshold
%             line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%                 reshape(nonzero_prob(i_setting,i_sim,num_threshold+k,i_cell,1:num_batch) ,[],1),...
%                 'Color',color_list(num_threshold+k,1,:),'LineStyle',line_styles{i_cell},'LineWidth',1)
%            end
           end
    end
    ylim([0 1]);
    xlabel('Number of trials','FontSize',15);
    ylabel('Non-zero probability','FontSize',25);
    title( strcat('Single batch \gamma_1 =', num2str(gamma_truth(1)),', \gamma_2=0.'),'FontSize',15); 
    hold off;
    
    figure(7)
    hold on;
    for i_cell = 1:n_cell
%            for k=1:num_threshold
%      
%             line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%                 reshape(mean_nonzero_prob(i_setting,1,k,i_cell,:) ,[],1),...
%                 'Color',color_list(k,1,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
%            end
%            for k=1:num_threshold
%             line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%                 reshape(mean_nonzero_prob(i_setting,1,num_threshold+k,i_cell,1:num_batch) ,[],1),...
%                 'Color',color_list(num_threshold+k,1,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
%            end
           for i_sim= 1:num_seed
%            for k=1:num_threshold
%             line((1: (trial_max/(K *n_cell)))*K*n_cell,...
%                 reshape(nonzero_prob(i_setting,i_sim,k,i_cell,:) ,[],1),...
%                 'Color',color_list(k,1,:),'LineStyle',line_styles{i_cell},'LineWidth',1)
%            end
           for k=1:num_threshold
            line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                reshape(nonzero_prob(i_setting,i_sim,num_threshold+k,i_cell,1:num_batch) ,[],1),...
                'Color',color_list(num_threshold+k,1,:),'LineStyle',line_styles{i_cell},'LineWidth',1)
           end
           end
    end
    ylim([0 1]);
    xlabel('Number of trials','FontSize',15);
    ylabel('Non-zero probability','FontSize',25);
    title( strcat('Full history \gamma_1 =', num2str(gamma_truth(1)),', \gamma_2=0.'),'FontSize',15); 
    hold off;
    
    saveas(1,strcat('./Figures/Oct18/gammaSim',num2str(i_setting),'.png'));
%     saveas(2,strcat('./Figures/Oct17/gainSim',num2str(i_setting),'.png'));
%     saveas(3,strcat('./Figures/Oct17/gamma_varianceSim',num2str(i_setting),'.png'));
%     saveas(4,strcat('./Figures/Oct17/gain_varianceSim',num2str(i_setting),'.png'));
%     saveas(5,strcat('./Figures/Oct11/TimeSim',num2str(i_setting),'.png'));
    saveas(6,strcat('./Figures/Oct18/ProbMass',num2str(i_setting),'.png'));
    saveas(7,strcat('./Figures/Oct18/ProbMass',num2str(i_setting),'Full.png'));
    
   
    close all;
end
