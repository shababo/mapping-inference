   num_trials_batch_grid   = [20   20  20  20  20  40  100 20  20  20  20  20  20]; 

for i_setting = 1:length(num_trials_batch_grid )
    N=4000/num_trials_batch_grid(i_setting);
    t_seq = [1 2 3 4 5*(1: N/5)];
    for design = 0:1
        for arg1 = 1:50
            
            flnm = strcat('./Data/Setting',num2str(i_setting),'Design', num2str(design));
            load(strcat(flnm, 'Rep',num2str(arg1),'_Truth.mat'));
            
            load(strcat(flnm, 'Rep',num2str(arg1),'_Crude.mat'));
            load(strcat(flnm, 'Rep',num2str(arg1),'_Gibbs.mat'));
            load(strcat(flnm, 'Rep',num2str(arg1),'_Sparse.mat')); %,'output_working');
            load(strcat(flnm, 'Rep',num2str(arg1),'_EM.mat')); %,'output_EM');
            
            % Read the estimates from the working model
            n_cell_local=size(Z,1);
            
            crude= cell(N,1);
            %threshold =  quantile([mpp_new.amplitudes], (1/num_threshold)*[0:num_threshold]);
            for i_batch = 1:N
                overall_connectivity = zeros(n_cell_local,1);
                overall_mark = zeros(n_cell_local,1);
                for j = 1:num_threshold
                    overall_connectivity = max(overall_connectivity, ...
                        output_crude{i_batch}(j).alpha(2:end));
                    for i_cell = 1:n_cell_local
                        if overall_connectivity(i_cell) == output_crude{i_batch}(j).alpha(i_cell+1)
                            overall_mark(i_cell)= (output_crude{i_batch}(num_threshold).threshold(j)+output_crude{i_batch}(num_threshold).threshold(j+1))/2;
                        end
                    end
                end
                crude{i_batch} = struct;
                crude{i_batch}.mean_gamma = overall_connectivity;
                crude{i_batch}.mean_mu = overall_mark;
                crude{i_batch}.comptime=time_record;
            end
            
            EM= cell(length(t_seq),1);
            for ind_t = 1:length(t_seq)
                EM{ind_t}.mean_gamma = output_EM{ind_t}.gamma_samples;
                EM{ind_t}.mean_mu = output_EM{ind_t}.mu_samples;
                EM{ind_t}.mean_sigma= output_EM{ind_t}.sigma_samples;
                EM{ind_t}.comptime=output_EM{ind_t}.delta_t ;
            end
            sparse= cell(length(t_seq),1);
            for ind_t = 1:length(t_seq)
                 sparse{ind_t}.mean_gamma = output_sparse{ind_t}.gamma_samples;
                 sparse{ind_t}.mean_mu = output_sparse{ind_t}.mu_samples;
                 sparse{ind_t}.mean_sigma= output_sparse{ind_t}.sigma_samples;
                 sparse{ind_t}.comptime=output_sparse{ind_t}.delta_t ;
            end
            
            %         Read the estimates from the Gibbs sampler (with soft assignments)
            Gibbs= cell(length(t_seq),1);
            for ind_t = 1:length(t_seq)
                Gibbs{ind_t}.mean_gamma =  mean(output_Gibbs{i_batch}.gamma_samples, 1);
                Gibbs{ind_t}.mean_mu =  mean(output_Gibbs{i_batch}.mu_samples, 1);
                Gibbs{ind_t}.mean_sigma =  mean(output_Gibbs{i_batch}.sigma_samples, 1);
                Gibbs{ind_t}.comptime = t_delta;
            end
            
            local_connected =local_amplitudes>0;
            % Output the relevant statistics:
            output= struct;
            n_cell_local= size(local_connected,1);
            dis_measure = zeros(n_cell_local,1);
            for i = 1:n_cell_local
                total_dist = ones(sum(local_connected),1)*Z(i,1:2)-Z(local_connected,1:2);
                sq_disc = total_dist.^2;
                dis_measure(i) = min(sqrt(sum(sq_disc,2)));
            end
            true_gamma =local_gamma;
            true_mu = local_amplitudes(local_gamma>0);
            
            output.dis_measure = dis_measure;
            output.local_amplitudes = local_amplitudes;
            
            output.local_connected = local_connected;
            output.true_gamma = true_gamma;
            output.true_mu = true_mu;
            output.crude=crude;
            
            output.EM= EM;
            output.sparse= sparse;
            output.Gibbs = Gibbs;
            
            
            outflnm=strcat('./Results/OutputDesign', num2str(design),'Setting',num2str(i_setting),...
                'Rep', num2str(arg1),'.mat');
            save(outflnm,'output');
        end
    end
end
%
%% Debugging:
 n_cell_local = size(output.true_gamma,1);
 jittered_gamma = output.true_gamma + normrnd(0,0.1,[n_cell_local 1]);
t_seq =[1 2 3 4 5*(1: N/5)];
   

 for ind_t = [1 2 3 10 12]

     figure(ind_t)
%      plot(jittered_gamma,output.EM{ind_t}.mean_gamma,'.','col',[0,1,0,1], 'markersize', 10);
%      hold on;
%      plot(jittered_gamma,output.working{ind_t}.mean_gamma,'.','col',[0,0,1,1], 'markersize', 10);
%      hold on;

     plot(jittered_gamma,output.crude{t_seq(ind_t)}.mean_gamma,'.','col',[1,0,0,1], 'markersize', 10);
     hold on;

     line([0 1],[0 1]);

     xlim([-0.1,1.1]);
     ylim([-0.1,1.1]);
 end
%
%
%%
%
% figure(11)
% temp = scatter(local_locations(:,1),...
%     -local_locations(:,2),...
%     (local_amplitudes+0.1)*35);
% set(temp,'MarkerFaceColor','k');
%     alpha(temp,0.3);
% hold on;
%     failed = scatter(local_locations(137,1),...
%     -local_locations(137,2),...
%     (local_amplitudes(137))*35);
% set(failed,'MarkerFaceColor','b');
%     alpha(failed,0.3);
%
%
% set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
% hold on;
% xlim([20,460]);
% ylim([-900,-400]);
%
% loc_trials = reshape([locations_trials], [(size(locations_trials,1))*4,1]);
% loc_trials_counts  = tabulate(loc_trials);
% loc_trials_counts_nonzero = loc_trials_counts(loc_trials_counts(:,2)>0,1:2);
%
% stimulate = scatter(Z_dense(loc_trials_counts_nonzero(:,1),1),...
%     -Z_dense(loc_trials_counts_nonzero(:,1),2),...
%     2*loc_trials_counts_nonzero(:,2), 'filled','d');
% set(stimulate,'MarkerFaceColor','r');
% alpha(stimulate,0.4);
% hold on;
