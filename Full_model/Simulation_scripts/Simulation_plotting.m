addpath(genpath('../../../mapping-inference'));
%% Specify the setting in this simulation
num_seed=20;
num_sim=6;
tpr_disconnected=zeros(num_sim,num_seed,100);
tpr_connected=zeros(num_sim,num_seed,100);
total_trial_learning=zeros(num_sim,num_seed,100);
total_trial_killing=zeros(num_sim,num_seed,100);
gamma_mean_path=zeros(num_sim,num_seed,100);
gamma_sd_path=zeros(num_sim,num_seed,100);
gamma_sd_path=zeros(num_sim,num_seed,100);

gain_mean_path=zeros(num_sim,num_seed,100);
gain_sd_path=zeros(num_sim,num_seed,100);
gain_sd_path=zeros(num_sim,num_seed,100);
iteration_counts=zeros(num_sim,num_seed);
gain_type=2;
good_prior=1;
for i_sim = 1:num_sim
    
    i_sim_index =i_sim+9;
    for i_seed = 1:num_seed
        % i_sim=0;i_seed=1;
        
        rng(i_seed,'twister');
        
        % Generate cellular parameters
        run('./Simulation_parameters.m')
        %% Load data
        load(strcat('./matfiles/Sep25/','Sim', num2str(i_sim_index),'Seed',num2str(i_seed),'.mat'))
        %'variational_params_path','gamma_path','var_gamma_path',...
        %'mpp_connected', 'trials_locations_connected','trials_powers_connected',...
        %'mpp_disconnected', 'trials_locations_disconnected','trials_powers_disconnected',...
        %'mpp_undefined', 'trials_locations_undefined','trials_powers_undefined',...
        %'undefined_cells', 'potentially_disconnected_cells', 'potentially_connected_cells',...
        %'dead_cells', 'alive_cells')
        %% Summarize results
        
        % Number of trials per group:
        final_iter = length(alive_cells);
        % Summarize the experiment:
        cell_assignments = zeros(final_iter,n_cell_this_plane);
        group_assignments = zeros(final_iter,2,5);
        cell_trial_numbers =zeros(final_iter,n_cell_this_plane); % number of trials on each cell
        cell_gamma_mean = zeros(final_iter,n_cell_this_plane);
        cell_gamma_variance = zeros(final_iter,n_cell_this_plane);
        cell_gain_mean = zeros(final_iter,n_cell_this_plane);
        cell_gain_variance = zeros(final_iter,n_cell_this_plane);
        trial_numbers = zeros(final_iter,3);
        connected_ind=find(gamma_truth(cell_group_list{this_plane})>0);
        disconnected_ind=find(gamma_truth(cell_group_list{this_plane})==0);
        total_trials_killing = zeros(final_iter,1);
        total_trials_learning = zeros(final_iter,1);
        
        for iter=2:final_iter
            % count the assignments
            cell_assignments(iter,:)= 5*alive_cells{iter}+4*potentially_connected_cells{iter}+...
                +undefined_cells{iter}*3+potentially_disconnected_cells{iter}*2+1*dead_cells{iter};
            for j = 1:5
                group_assignments(iter,2,j)= sum(cell_assignments(iter,connected_ind)==j);
                group_assignments(iter,1,j)= sum(cell_assignments(iter,disconnected_ind)==j);
            end
            
            % count the number of trials
            
            % map the stim location to cells
            
            %     for j=1:n_cell_this_plane
            %     cell_trial_numbers(iter,j)= sum(loc_to_cell_nuclei([trials_locations_connected{iter-1}])==j)+...
            %         sum(sum(trials_locations_undefined{iter-1}==j))+sum(sum(trials_locations_disconnected{iter-1}==j));
            %     end
            trial_numbers(iter,3)=size(trials_locations_connected{iter-1},1);
            trial_numbers(iter,2)=size(trials_locations_undefined{iter-1},1);
            trial_numbers(iter,1)=size(trials_locations_disconnected{iter-1},1);
            
            % record the gamma and gain
            for j=1:n_cell_this_plane
                [cell_gain_mean(iter,j), cell_gain_variance(iter,j)]=calculate_posterior_mean(...
                    variational_params_path.alpha_gain(j,iter),variational_params_path.beta_gain(j,iter),gain_bound.low,gain_bound.up);
                [cell_gamma_mean(iter,j), cell_gamma_variance(iter,j)]=calculate_posterior_mean(...
                    variational_params_path.alpha(j,iter),variational_params_path.beta(j,iter),0,1);
            end
        end
        
        
        % Record the useful info: 
        tpr_disconnected(i_sim+1,i_seed,1:final_iter)=reshape(group_assignments(:,1,1),[1 1 final_iter])/length(disconnected_ind);
        tpr_connected(i_sim+1,i_seed,1:final_iter)=reshape(group_assignments(:,2,5),[1 1 final_iter])/length(connected_ind);
        total_trial_killing(i_sim+1,i_seed,1:final_iter)= cumsum(sum(trial_numbers(:,1:2),2));
        total_trial_learning(i_sim+1,i_seed,1:final_iter)= cumsum(trial_numbers(:,3));
        gamma_mean_path(i_sim+1,i_seed,1:final_iter)=mean(abs(cell_gamma_mean(:,connected_ind)-gamma_truth(target_cell_list(1).primary(connected_ind))'),2);
        gamma_sd_path(i_sim+1,i_seed,1:final_iter)=sqrt(mean(cell_gamma_variance(:,connected_ind),2));
        gamma_sd_path(i_sim+1,i_seed,1)=1;

        gain_mean_path(i_sim+1,i_seed,1:final_iter)=mean(abs(cell_gain_mean(:,connected_ind)-gain_truth(target_cell_list(1).primary(connected_ind))')./gain_truth(target_cell_list(1).primary(connected_ind))',2);
        gain_sd_path(i_sim+1,i_seed,1:final_iter)=sqrt(mean(cell_gain_variance(:,connected_ind),2));
        gain_sd_path(i_sim+1,i_seed,1)=1;
        
        iteration_counts(i_sim+1,i_seed)=final_iter;
            
    end
end

%% Average results:
n_grid=3000;
trial_num_grid = 1:n_grid;

tpr_disconnected_mean=zeros(num_sim+1,n_grid);
tpr_connected_mean=zeros(num_sim+1,n_grid);
gamma_mean_path_mean=zeros(num_sim+1,n_grid);
gamma_sd_path_mean=zeros(num_sim+1,n_grid);
gamma_sd_path_mean=zeros(num_sim+1,n_grid);
gain_mean_path_mean=zeros(num_sim+1,n_grid);
gain_sd_path_mean=zeros(num_sim+1,n_grid);

for i_sim = 1:num_sim
    for i_seed = 1:num_seed
        [unique_trial,i_loc,~]=unique(reshape(total_trial_killing(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)), [],1));
        voi=reshape(tpr_disconnected(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)),[],1);
        voi=voi(i_loc);
        temp_mean = interp1(unique_trial,...
            voi,...
            trial_num_grid);
        temp_mean(isnan(temp_mean))=temp_mean(min(find(isnan(temp_mean))-1));
        tpr_disconnected_mean(i_sim+1,:)=tpr_disconnected_mean(i_sim+1,:)+temp_mean/num_seed;
        
        
        [unique_trial,i_loc,~]=unique(reshape(total_trial_learning(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)), [],1));
        voi= reshape(tpr_connected(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)),[],1);
        voi=voi(i_loc);
        temp_mean = interp1(unique_trial,...
            voi,...
            trial_num_grid);
        temp_mean(isnan(temp_mean))=temp_mean(min(find(isnan(temp_mean))-1));
        tpr_connected_mean(i_sim+1,:)=tpr_connected_mean(i_sim+1,:)+temp_mean/num_seed;
        
        
        [unique_trial,i_loc,~]=unique(reshape(total_trial_learning(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)), [],1));
        voi= reshape(gamma_mean_path(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)),[],1);
        voi=voi(i_loc);
        temp_mean = interp1(unique_trial,...
            voi,...
            trial_num_grid);
        temp_mean(isnan(temp_mean))=temp_mean(min(find(isnan(temp_mean))-1)) ;
        gamma_mean_path_mean(i_sim+1,:)=gamma_mean_path_mean(i_sim+1,:)+temp_mean/num_seed;
        
        [unique_trial,i_loc,~]=unique(reshape(total_trial_learning(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)), [],1));
        voi= reshape(gamma_sd_path(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)),[],1);
        voi=voi(i_loc);
        temp_mean = interp1(unique_trial,...
            voi,...
           trial_num_grid);
        temp_mean(isnan(temp_mean))=temp_mean(min(find(isnan(temp_mean))-1));
        gamma_sd_path_mean(i_sim+1,:)=gamma_sd_path_mean(i_sim+1,:)+temp_mean/num_seed;
        
        [unique_trial,i_loc,~]=unique(reshape(total_trial_learning(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)), [],1));
        voi= reshape(gain_mean_path(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)),[],1);
        voi=voi(i_loc);
         temp_mean = interp1(unique_trial,...
            voi,...
           trial_num_grid);
        temp_mean(isnan(temp_mean))=temp_mean(min(find(isnan(temp_mean))-1));
        gain_mean_path_mean(i_sim+1,:)=gain_mean_path_mean(i_sim+1,:)+temp_mean/num_seed;
        
        [unique_trial,i_loc,~]=unique(reshape(total_trial_learning(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)), [],1));
        voi= reshape(gain_sd_path(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)),[],1);
        voi=voi(i_loc);
         temp_mean = interp1(unique_trial,...
            voi,...
            trial_num_grid);
        temp_mean(isnan(temp_mean))=temp_mean(min(find(isnan(temp_mean))-1));
        gain_sd_path_mean(i_sim+1,:)=gain_sd_path_mean(i_sim+1,:)+temp_mean/num_seed;
    end
end
%% Plotting
% Sample size v.s. TPR (for disconnected cells)
color_list= [[1 0 0 0.1]; [1 0 0 0.1]; [0 0 0 0.1]; [0 0 1 0.1]; [0 1 0 0.1]];
color_list_solid= [[1 0 0 1]; [1 0 0 1]; [0 0 0 1]; [0 0 1 1]; [0 1 0 1]];
i_sim=5;
% for i_sim = 1:4
figure(1)
hold on;
for i_seed = 1:num_seed
line(reshape(total_trial_killing(i_sim+2,i_seed,1:iteration_counts(i_sim+1,i_seed)), [],1),...
    reshape(tpr_disconnected(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)),[],1),...
    'Color',color_list(1,:),...
    'LineStyle','-','LineWidth',1)
line(reshape(total_trial_killing(i_sim+1,i_seed,1:iteration_counts(1,i_seed)),[],1),...
    reshape(tpr_disconnected(1,i_seed,1:iteration_counts(1,i_seed)),[],1),...
    'Color',color_list(1,:),...
    'LineStyle',':','LineWidth',1)

end
        
line(trial_num_grid,...
    tpr_disconnected_mean(i_sim+2,:),...
    'Color',color_list_solid(1,:),...
    'LineStyle','-','LineWidth',3)
line(trial_num_grid,...
    tpr_disconnected_mean(i_sim+1,:),...
    'Color',color_list_solid(1,:),...
    'LineStyle',':','LineWidth',3)
xlim([0 600])
xlabel('Number of trials on pot. disconnected and undefined cells','FontSize',15);
ylabel('True negative rate','FontSize',25);
hold off;


figure(2)
hold on;for i_seed = 1:num_seed
line(reshape(total_trial_learning(i_sim+2,i_seed,1:iteration_counts(i_sim+1,i_seed)), [],1),...
    reshape(tpr_connected(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)),[],1),...
    'Color',color_list(4,:),...
    'LineStyle','-','LineWidth',1)

line(reshape(total_trial_learning(i_sim+1,i_seed,1:iteration_counts(1,i_seed)),[],1),...
    reshape(tpr_connected(1,i_seed,1:iteration_counts(1,i_seed)),[],1),...
    'Color',color_list(4,:),...
    'LineStyle',':','LineWidth',1)
end
line(trial_num_grid,...
    tpr_connected_mean(i_sim+2,:),...
    'Color',color_list_solid(4,:),...
    'LineStyle','-','LineWidth',3)
line(trial_num_grid,...
    tpr_connected_mean(i_sim+1,:),...
    'Color',color_list_solid(4,:),...
    'LineStyle',':','LineWidth',3)

xlabel('Number of trials on pot. connected cells','FontSize',15);
ylabel('True positive rate','FontSize',25);
hold off;


figure(3)
hold on;
for i_seed = 1:num_seed
line(reshape(total_trial_learning(i_sim+2,i_seed,1:iteration_counts(i_sim+1,i_seed)), [],1),...
    reshape(gamma_mean_path(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)),[],1),...
    'Color',color_list(4,:),...
    'LineStyle','-','LineWidth',1)
line(reshape(total_trial_learning(i_sim+1,i_seed,1:iteration_counts(1,i_seed)),[],1),...
    reshape(gamma_mean_path(1,i_seed,1:iteration_counts(1,i_seed)),[],1),...
    'Color',color_list(4,:),...
    'LineStyle',':','LineWidth',1)
end

line(trial_num_grid,...
    gamma_mean_path_mean(i_sim+2,:),...
    'Color',color_list_solid(4,:),...
    'LineStyle','-','LineWidth',3)
line(trial_num_grid,...
    gamma_mean_path_mean(i_sim+1,:),...
    'Color',color_list_solid(4,:),...
    'LineStyle',':','LineWidth',3)
xlabel('Number of trials on potentially connected cells','FontSize',15);
ylabel('|E[\gamma|data]-\gamma|','FontSize',25);
hold off;


figure(4)
hold on;
for i_seed = 1:num_seed
line(reshape(total_trial_learning(i_sim+2,i_seed,1:iteration_counts(i_sim+1,i_seed)), [],1),...
    reshape(gamma_sd_path(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)),[],1),...
    'Color',color_list(5,:),...
    'LineStyle','-','LineWidth',1)
line(reshape(total_trial_learning(i_sim+1,i_seed,1:iteration_counts(1,i_seed)),[],1),...
    reshape(gamma_sd_path(1,i_seed,1:iteration_counts(1,i_seed)),[],1),...
    'Color',color_list(5,:),...
    'LineStyle',':','LineWidth',1)
end
line(trial_num_grid,...
    gamma_sd_path_mean(i_sim+2,:),...
    'Color',color_list_solid(5,:),...
    'LineStyle','-','LineWidth',3)
line(trial_num_grid,...
    gamma_sd_path_mean(i_sim+1,:),...
    'Color',color_list_solid(5,:),...
    'LineStyle',':','LineWidth',3)

xlabel('Number of trials on potentially connected cells','FontSize',15);
ylabel('sd[\gamma|data]','FontSize',25);
hold off;

figure(5)
hold on;
for i_seed = 1:num_seed
line(reshape(total_trial_learning(i_sim+2,i_seed,1:iteration_counts(i_sim+1,i_seed)), [],1),...
    reshape(gain_mean_path(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)),[],1),...
    'Color',color_list(4,:),...
    'LineStyle','-','LineWidth',1)
line(reshape(total_trial_learning(i_sim+1,i_seed,1:iteration_counts(1,i_seed)),[],1),...
    reshape(gain_mean_path(1,i_seed,1:iteration_counts(1,i_seed)),[],1),...
    'Color',color_list(4,:),...
    'LineStyle',':','LineWidth',1)
end

line(trial_num_grid,...
    gain_mean_path_mean(i_sim+2,:),...
    'Color',color_list_solid(4,:),...
    'LineStyle','-','LineWidth',3)
line(trial_num_grid,...
    gain_mean_path_mean(i_sim+1,:),...
    'Color',color_list_solid(4,:),...
    'LineStyle',':','LineWidth',3)
hold off;

figure(6)
hold on;
for i_seed = 1:num_seed
line(reshape(total_trial_learning(i_sim+2,i_seed,1:iteration_counts(i_sim+1,i_seed)), [],1),...
    reshape(gain_sd_path(i_sim+1,i_seed,1:iteration_counts(i_sim+1,i_seed)),[],1),...
    'Color',color_list(5,:),...
    'LineStyle','-','LineWidth',1)
line(reshape(total_trial_learning(i_sim+1,i_seed,1:iteration_counts(1,i_seed)),[],1),...
    reshape(gain_sd_path(1,i_seed,1:iteration_counts(1,i_seed)),[],1),...
    'Color',color_list(5,:),...
    'LineStyle',':','LineWidth',1)
end

line(trial_num_grid,...
    gain_sd_path_mean(i_sim+2,:),...
    'Color',color_list_solid(5,:),...
    'LineStyle','-','LineWidth',3)
line(trial_num_grid,...
    gain_sd_path_mean(i_sim+1,:),...
    'Color',color_list_solid(5,:),...
    'LineStyle',':','LineWidth',3)
hold off;

saveas(1,strcat('./Figures/Sep25/Sim',num2str(i_sim),'tnr_trials','.png'));
saveas(2,strcat('./Figures/Sep25/Sim',num2str(i_sim),'tpr_trials','.png'));
saveas(3,strcat('./Figures/Sep25/Sim',num2str(i_sim),'gamma_mean_trials','.png'));
saveas(4,strcat('./Figures/Sep25/Sim',num2str(i_sim),'gamma_sd_trials','.png'));
  close all
% end

%% New methods 
color_list= [[1 0 0 0.1]; [0 1 0 0.1]; [0 0 1 0.1]; [0 0 1 0.1]; [0 1 0 0.1]];
color_list_solid= [[1 0 0 1]; [0 1 0 1]; [0 0 1 1]; [0 0 1 1]; [0 1 0 1]];

for i_setting = 1:2
figure(1)
hold on;
for i_method = 1:3
for i_seed = 1:num_seed
line(reshape(total_trial_killing( 3*(i_setting-1)+i_method+1,i_seed,1:iteration_counts( 3*(i_setting-1)+i_method+1,i_seed)), [],1),...
    reshape(tpr_disconnected( 3*(i_setting-1)+i_method+1,i_seed,1:iteration_counts( 3*(i_setting-1)+i_method+1,i_seed)),[],1),...
    'Color',color_list(i_method,:),...
    'LineStyle','-','LineWidth',1)
end
line(trial_num_grid,...
    tpr_disconnected_mean(3*(i_setting-1)+i_method+1,:),...
    'Color',color_list_solid(i_method,:),...
    'LineStyle','-','LineWidth',3)
end
xlim([0 1000])
xlabel('Number of trials on pot. disconnected and undefined cells','FontSize',15);
ylabel('True negative rate','FontSize',25);
hold off;


figure(2)
hold on;

for i_method = 1:3
for i_seed = 1:num_seed
line(reshape(total_trial_learning( 3*(i_setting-1)+i_method+1,i_seed,1:iteration_counts( 3*(i_setting-1)+i_method+1,i_seed)), [],1),...
    reshape(tpr_connected( 3*(i_setting-1)+i_method+1,i_seed,1:iteration_counts( 3*(i_setting-1)+i_method+1,i_seed)),[],1),...
    'Color',color_list(i_method,:),...
    'LineStyle','-','LineWidth',1)
end
line(trial_num_grid,...
    tpr_connected_mean(3*(i_setting-1)+i_method+1,:),...
    'Color',color_list_solid(i_method,:),...
    'LineStyle','-','LineWidth',3)
end

xlabel('Number of trials on pot. connected cells','FontSize',15);
ylabel('True positive rate','FontSize',25);
hold off;


figure(3)
hold on;

for i_method = 1:3
    for i_seed = 1:num_seed
    line(reshape(total_trial_learning( 3*(i_setting-1)+i_method+1,i_seed,1:iteration_counts( 3*(i_setting-1)+i_method+1,i_seed)), [],1),...
        reshape(gamma_mean_path(3*(i_setting-1)+i_method+1,1:iteration_counts( 3*(i_setting-1)+i_method+1,i_seed)),[],1),...
        'Color',color_list(i_method,:),...
        'LineStyle','-','LineWidth',1)
    end
line(trial_num_grid,...
    gamma_mean_path_mean(3*(i_setting-1)+i_method+1,:),...
    'Color',color_list_solid(i_method,:),...
    'LineStyle','-','LineWidth',3)
end
xlabel('Number of trials on potentially connected cells','FontSize',15);
ylabel('|E[\gamma|data]-\gamma|','FontSize',25);
hold off;



figure(4)
hold on;

for i_method = 1:3
for i_seed = 1:num_seed
line(reshape(total_trial_learning( 3*(i_setting-1)+i_method+1,i_seed,1:iteration_counts( 3*(i_setting-1)+i_method+1,i_seed)), [],1),...
    reshape( gamma_sd_path( 3*(i_setting-1)+i_method+1,i_seed,1:iteration_counts( 3*(i_setting-1)+i_method+1,i_seed)),[],1),...
    'Color',color_list(i_method,:),...
    'LineStyle','-','LineWidth',1)
end
line(trial_num_grid,...
    gamma_sd_path_mean(3*(i_setting-1)+i_method+1,:),...
    'Color',color_list_solid(i_method,:),...
    'LineStyle','-','LineWidth',3)
end
xlabel('Number of trials on potentially connected cells','FontSize',15);
ylabel('sd[\gamma|data]','FontSize',25);
hold off;

saveas(1,strcat('./Figures/Sep25/Setting',num2str(i_setting),'tnr_trials','.png'));
saveas(2,strcat('./Figures/Sep25/Setting',num2str(i_setting),'tpr_trials','.png'));
saveas(3,strcat('./Figures/Sep25/Setting',num2str(i_setting),'gamma_mean_trials','.png'));
saveas(4,strcat('./Figures/Sep25/Setting',num2str(i_setting),'gamma_sd_trials','.png'));
  close all
end

%%
figure(5)
for i_method = 1:3
    for i_seed = 1:num_seed
    line(reshape(total_trial_learning( 3*(i_setting-1)+i_method+1,i_seed,2:iteration_counts( 3*(i_setting-1)+i_method+1,i_seed)), [],1),...
        reshape(gain_mean_path(3*(i_setting-1)+i_method+1,2:iteration_counts( 3*(i_setting-1)+i_method+1,i_seed)),[],1),...
        'Color',color_list(i_method,:),...
        'LineStyle','-','LineWidth',1)
    end
    line(trial_num_grid,...
        gain_mean_path_mean(3*(i_setting-1)+i_method+1,:),...
        'Color',color_list_solid(i_method,:),...
        'LineStyle','-','LineWidth',3)
end
ylim([0, 1]);
xlabel('Number of trials on potentially connected cells','FontSize',15);
ylabel('|E[\phi|data]-\phi|','FontSize',25);
hold off;

figure(6)
hold on;
for i_method = 1:3
    for i_seed = 1:num_seed
    line(reshape(total_trial_learning( 3*(i_setting-1)+i_method+1,i_seed,2:iteration_counts( 3*(i_setting-1)+i_method+1,i_seed)), [],1),...
        reshape(gain_sd_path(3*(i_setting-1)+i_method+1,2:iteration_counts( 3*(i_setting-1)+i_method+1,i_seed)),[],1),...
        'Color',color_list(i_method,:),...
        'LineStyle','-','LineWidth',1)
    end
    line(trial_num_grid,...
        gain_sd_path_mean(3*(i_setting-1)+i_method+1,:),...
        'Color',color_list_solid(i_method,:),...
        'LineStyle','-','LineWidth',3)
end
ylim([0, 0.01]);
xlabel('Number of trials on potentially connected cells','FontSize',15);
% ylabel('|E[\phi|data]-\phi|','FontSize',25);
hold off;
