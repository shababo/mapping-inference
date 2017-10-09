addpath(genpath('../../../mapping-inference'));
%% Parameters:
d=10; % distance between the two cells
%% Cell locations
cell_locations=zeros(3,3);
cell_locations(2,1)=d;
cell_locations(3,1)=d/2;
cell_locations(3,2)=sqrt(3)*d/2;
n_cell=size(cell_locations,1);
%% Cell parameters
num_sim=3;
num_seed=80;

for i_setting = 1:num_sim
    background_rate=1e-4;
    v_th_known=15*ones([n_cell,1]);v_reset_known=-1e4*ones([n_cell,1]);
    g_truth = 0.02*ones([n_cell,1]);
    funcs.invlink = @invlink_sig;%@(resp) log(1 + exp(resp));%@(x) exp(x);
    if i_setting == 1
        gamma_truth = [0.6;0.7;0.8];gain_truth=[0.015; 0.02; 0.025];
    elseif i_setting == 2
        gamma_truth = [0.6;0.7;0];gain_truth=[0.015; 0.02; 0.025];
    elseif i_setting == 3
        gamma_truth = [0.6;0;0];gain_truth=[0.015; 0.02; 0.025];
    end
    
    
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
        for i_cell = 1:n_cell
            loc_to_cell( (2*num_per_grid+1)*(i_cell-1)+ (1:(2*num_per_grid+1) ))=i_cell;
        end
        power_selected=zeros(n_cell*(2*num_per_grid+1),1);
        power_sd=zeros(n_cell*(2*num_per_grid+1),1);
        
        %% Parameters in the design stage
        % Design parameters
        n_spots_per_trial = 1;
        % Need to run sims to check how these parameters affect the results
        % Prior distribution
        prior_pi0=0.8;
        
        % Initialize tuning parameters in the VI
        C_threshold = 0.01;maxit=1000;
        S=200;epsilon=0.01;eta_logit=0;
        background_rt=background_rate*time_max; % raw probability of firing within a trial
        eta_beta=0.05;
        
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
        
        var_gain_prior=1;
        gain_bias=0;
        K=5;
        n_replicates=1;
        cell_list= 1:n_cell;
        trial_max=300;
        iter=1;
        n_MC_samples=100;
        epislon=10;
        connected_threshold=0.4;
        %%
        run('Twoneurons_naive.m');
        n_events_naive
        
        %%
        connected_threshold=1;
        run('Twoneurons_optimal.m');
        n_events_optimal
        gamma_path_pure=gamma_path;
        gain_path_pure=gain_path;
        gamma_25_path_pure=gamma_25_path;
        gamma_75_path_pure=gamma_75_path;
        gain_25_path_pure=gain_25_path;
        gain_75_path_pure=gain_75_path;
        %%
        connected_threshold=0;
        run('Twoneurons_optimal.m');
        n_events_optimal
        gamma_path_jit=gamma_path;
        gain_path_jit=gain_path;
        gamma_25_path_jit=gamma_25_path;
        gamma_75_path_jit=gamma_75_path;
        gain_25_path_jit=gain_25_path;
        gain_75_path_jit=gain_75_path;
        
        %%
        save(strcat('./matfiles/Oct04/','Sim', num2str(i_setting),'Seed',num2str(i_seed),'.mat'),...
            'gamma_path_naive','gamma_25_path_naive','gamma_75_path_naive',...
            'gain_path_naive','gain_25_path_naive','gain_75_path_naive',...
            'mpp_naive','trials_locations_naive','trials_powers_naive',...
            'gamma_path_optimal','gamma_25_path_optimal','gamma_75_path_optimal',...
            'gain_path_optimal','gain_25_path_optimal','gain_75_path_optimal',...
            'gamma_path_pure','gamma_25_path_pure','gamma_75_path_pure',...
            'gain_path_pure','gain_25_path_pure','gain_75_path_pure',...
            'gamma_path_jit','gamma_25_path_jit','gamma_75_path_jit',...
            'gain_path_jit','gain_25_path_jit','gain_75_path_jit',...
            'mpp_optimal','trials_locations_optimal','trials_powers_optimal',...
            'K','trial_max','n_cell')
    end
end
%%
num_sim=3;
num_seed=80;
n_cell=3;
K=5;
trial_max=300;
gamma_all_path_naive=zeros(num_sim,num_seed,3,n_cell,(trial_max/(K *n_cell)));
gain_all_path_naive=zeros(num_sim,num_seed,3,n_cell,(trial_max/(K *n_cell)));
gamma_deviation_naive=zeros(num_sim,num_seed,n_cell,(trial_max/(K *n_cell)));
gain_deviation_naive=zeros(num_sim,num_seed,n_cell,(trial_max/(K *n_cell)));

gamma_all_path_optimal=zeros(num_sim,num_seed,3,n_cell,(trial_max/(K *n_cell)));
gain_all_path_optimal=zeros(num_sim,num_seed,3,n_cell,(trial_max/(K *n_cell)));
gamma_deviation_optimal=zeros(num_sim,num_seed,n_cell,(trial_max/(K *n_cell)));
gain_deviation_optimal=zeros(num_sim,num_seed,n_cell,(trial_max/(K *n_cell)));

gamma_all_path_pure=zeros(num_sim,num_seed,3,n_cell,(trial_max/(K *n_cell)));
gain_all_path_pure=zeros(num_sim,num_seed,3,n_cell,(trial_max/(K *n_cell)));
gamma_deviation_pure=zeros(num_sim,num_seed,n_cell,(trial_max/(K *n_cell)));
gain_deviation_pure=zeros(num_sim,num_seed,n_cell,(trial_max/(K *n_cell)));

gamma_all_path_jit=zeros(num_sim,num_seed,3,n_cell,(trial_max/(K *n_cell)));
gain_all_path_jit=zeros(num_sim,num_seed,3,n_cell,(trial_max/(K *n_cell)));
gamma_deviation_jit=zeros(num_sim,num_seed,n_cell,(trial_max/(K *n_cell)));
gain_deviation_jit=zeros(num_sim,num_seed,n_cell,(trial_max/(K *n_cell)));


%
for i_setting = 1:num_sim
    if i_setting == 1
        gamma_truth = [0.6;0.7;0.8];gain_truth=[0.015; 0.02; 0.025];
    elseif i_setting == 2
        gamma_truth = [0.6;0.7;0];gain_truth=[0.015; 0.02; 0.025];
    elseif i_setting == 3
        gamma_truth = [0.6;0;0];gain_truth=[0.015; 0.02; 0.025];
    end
    for i_seed = 1:num_seed
        % i_sim=0;i_seed=1;
        
        % Load data
        load(strcat('./matfiles/Oct04/','Sim', num2str(i_setting),'Seed',num2str(i_seed),'.mat'))
        gamma_all_path_naive(i_setting,i_seed,1,:,:)=gamma_path_naive(:,2:end);
        gain_all_path_naive(i_setting,i_seed,1,:,:)=gain_path_naive(:,2:end);
        gamma_deviation_naive(i_setting,i_seed,:,:)=(gamma_path_naive(:,2:end)-gamma_truth).^2;
        gain_deviation_naive(i_setting,i_seed,:,:)=(gain_path_naive(:,2:end)-gain_truth).^2;
        
        gamma_all_path_optimal(i_setting,i_seed,1,:,:)=gamma_path_optimal(:,2:end);
        gain_all_path_optimal(i_setting,i_seed,1,:,:)=gain_path_optimal(:,2:end);
        gamma_deviation_optimal(i_setting,i_seed,:,:)=(gamma_path_optimal(:,2:end)-gamma_truth).^2;
        gain_deviation_optimal(i_setting,i_seed,:,:)=(gain_path_optimal(:,2:end)-gain_truth).^2;
        
        
        
        gamma_all_path_pure(i_setting,i_seed,1,:,:)=gamma_path_pure(:,2:end);
        gain_all_path_pure(i_setting,i_seed,1,:,:)=gain_path_pure(:,2:end);
        gamma_deviation_pure(i_setting,i_seed,:,:)=(gamma_path_pure(:,2:end)-gamma_truth).^2;
        gain_deviation_pure(i_setting,i_seed,:,:)=(gain_path_pure(:,2:end)-gain_truth).^2;
        
        
        gamma_all_path_jit(i_setting,i_seed,1,:,:)=gamma_path_jit(:,2:end);
        gain_all_path_jit(i_setting,i_seed,1,:,:)=gain_path_jit(:,2:end);
        gamma_deviation_jit(i_setting,i_seed,:,:)=(gamma_path_jit(:,2:end)-gamma_truth).^2;
        gain_deviation_jit(i_setting,i_seed,:,:)=(gain_path_jit(:,2:end)-gain_truth).^2;
    end
end
%% Take average

mean_gamma_all_path_naive=mean(gamma_all_path_naive,2);
mean_gain_all_path_naive=mean(gain_all_path_naive,2);
mean_gain_deviation_naive=mean(gain_deviation_naive,2);
mean_gamma_deviation_naive=mean(gamma_deviation_naive,2);

mean_gamma_all_path_optimal=mean(gamma_all_path_optimal,2);
mean_gain_all_path_optimal=mean(gain_all_path_optimal,2);
mean_gain_deviation_optimal=mean(gain_deviation_optimal,2);
mean_gamma_deviation_optimal=mean(gamma_deviation_optimal,2);

mean_gamma_all_path_pure=mean(gamma_all_path_pure,2);
mean_gain_all_path_pure=mean(gain_all_path_pure,2);
mean_gain_deviation_pure=mean(gain_deviation_pure,2);
mean_gamma_deviation_pure=mean(gamma_deviation_pure,2);

mean_gamma_all_path_jit=mean(gamma_all_path_jit,2);
mean_gain_all_path_jit=mean(gain_all_path_jit,2);
mean_gain_deviation_jit=mean(gain_deviation_jit,2);
mean_gamma_deviation_jit=mean(gamma_deviation_jit,2);


%% Plot:
line_styles={'-' ':' '--'};
color_list_naive=[ [1 0 0 0.7];[1 0 0 0.3];[1 0 0 0.3]];
color_list_optimal=[ [0 0 0 0.7];[0 0 0 0.3];[0 0 0 0.3]];
color_list_pure=[ [0 0 1 0.7];[0 0 1 0.3];[0 0 1 0.3]];
color_list_jit=[ [1 0 1 0.7];[1 0 1 0.3];[1 0 1 0.3]];

for i_setting = 1:num_sim
    if i_setting == 1
        gamma_truth = [0.6;0.7;0.8];gain_truth=[0.015; 0.02; 0.025];
    elseif i_setting == 2
        gamma_truth = [0.6;0.7;0];gain_truth=[0.015; 0.02; 0.025];
    elseif i_setting == 3
        gamma_truth = [0.6;0;0];gain_truth=[0.015; 0.02; 0.025];
    end
    
    % i_sim=0;i_seed=1;
    
    figure(1)
    hold on;
    for i_cell = 1:n_cell
        for k=1:1
            line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                reshape(mean_gamma_all_path_naive(i_setting,1,k,i_cell,:) ,[],1),...
                'Color',color_list_naive(k,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
        end
    end
    for i_cell = 1:n_cell
        for k=1:1
            line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                reshape(mean_gamma_all_path_optimal(i_setting,1,k,i_cell,:) ,[],1),...
                'Color',color_list_optimal(k,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
        end
    end
    
    for i_cell = 1:n_cell
        for k=1:1
            line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                reshape(mean_gamma_all_path_pure(i_setting,1,k,i_cell,:) ,[],1),...
                'Color',color_list_pure(k,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
        end
    end
    
    for i_cell = 1:n_cell
        for k=1:1
            line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                reshape(mean_gamma_all_path_jit(i_setting,1,k,i_cell,:) ,[],1),...
                'Color',color_list_jit(k,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
        end
    end
    
    for i_cell = 1:n_cell
        line([K*n_cell trial_max], [gamma_truth(i_cell) gamma_truth(i_cell)],...
            'Color',[0 1 0 0.5],...
            'LineStyle',line_styles{i_cell},'LineWidth',3)
    end
    ylim([0 1]);
    xlabel('Number of trials','FontSize',15);
    ylabel('E[\gamma|data]','FontSize',25);
    hold off;
    %
    
    
    figure(2)
    
    hold on;
    for i_cell = 1:n_cell
        if gamma_truth(i_cell)>0
            
            for k=1:1
                line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                    reshape(mean_gain_all_path_naive(i_setting,1,k,i_cell,:) ,[],1),...
                    'Color',color_list_naive(k,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
                
                line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                    reshape(mean_gain_all_path_optimal(i_setting,1,k,i_cell,:) ,[],1),...
                    'Color',color_list_optimal(k,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
                
                line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                    reshape(mean_gain_all_path_pure(i_setting,1,k,i_cell,:) ,[],1),...
                    'Color',color_list_pure(k,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
                
                line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                    reshape(mean_gain_all_path_jit(i_setting,1,k,i_cell,:) ,[],1),...
                    'Color',color_list_jit(k,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
                
                line([K*n_cell trial_max], [gain_truth(i_cell) gain_truth(i_cell)],...
                    'Color',[0 1 0 0.5],...
                    'LineStyle',line_styles{i_cell},'LineWidth',3)
            end
        end
    end
    
    
    xlabel('Number of trials','FontSize',15);
    ylabel('E[\phi|data]','FontSize',25);
    hold off;
    
    figure(3)
    hold on;
    for i_cell = 1:n_cell
        for k=1:1
            line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                reshape(mean_gamma_deviation_naive(i_setting,1,i_cell,:) ,[],1),...
                'Color',color_list_naive(k,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
        end
    end
    for i_cell = 1:n_cell
        for k=1:1
            line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                reshape(mean_gamma_deviation_optimal(i_setting,1,i_cell,:) ,[],1),...
                'Color',color_list_optimal(k,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
        end
    end
    for i_cell = 1:n_cell
        for k=1:1
            line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                reshape(mean_gamma_deviation_pure(i_setting,1,i_cell,:) ,[],1),...
                'Color',color_list_pure(k,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
        end
    end
    
    for i_cell = 1:n_cell
        for k=1:1
            line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                reshape(mean_gamma_deviation_jit(i_setting,1,i_cell,:) ,[],1),...
                'Color',color_list_jit(k,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
        end
    end
    xlabel('Number of trials','FontSize',15);
    ylabel('||E[\gamma|data]-\gamma^*||','FontSize',25);
    hold off;
    
    
    figure(4)
    hold on;
    for i_cell = 1:n_cell
        if gamma_truth(i_cell)>0
            
            for k=1:1
                line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                    reshape(mean_gain_deviation_naive(i_setting,1,i_cell,:) ,[],1),...
                    'Color',color_list_naive(k,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
                line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                    reshape(mean_gain_deviation_optimal(i_setting,1,i_cell,:) ,[],1),...
                    'Color',color_list_optimal(k,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
                line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                    reshape(mean_gain_deviation_pure(i_setting,1,i_cell,:) ,[],1),...
                    'Color',color_list_pure(k,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
                line((1: (trial_max/(K *n_cell)))*K*n_cell,...
                    reshape(mean_gain_deviation_jit(i_setting,1,i_cell,:) ,[],1),...
                    'Color',color_list_jit(k,:),'LineStyle',line_styles{i_cell},'LineWidth',3)
            end
        end
    end
    xlabel('Number of trials','FontSize',15);
    ylabel('||E[\phi|data]-\phi^*||','FontSize',25);
    hold off;
    
    saveas(1,strcat('./Figures/Oct04/VertSim',num2str(i_setting),'gamma','.png'));
    saveas(2,strcat('./Figures/Oct04/VertSim',num2str(i_setting),'gain','.png'));
    saveas(3,strcat('./Figures/Oct04/VertSim',num2str(i_setting),'gamma_mse','.png'));
    saveas(4,strcat('./Figures/Oct04/VertSim',num2str(i_setting),'gain_mse','.png'));
    
    close all;
end
