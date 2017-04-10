temp_struct = struct;
temp_struct.crude=zeros(N,1);
temp_struct.Gibbs=zeros(length(t_seq),1);
temp_struct.EM=zeros(length(t_seq),1);

summary_stats = cell(2);
i_rep =  [1:50];
n_rep=length(i_rep);
for design = 0:1
    NRE=temp_struct;
    AUC=temp_struct;
    NRE_mu=temp_struct;
    
    for arg1 = i_rep
        flnm=strcat('./Results/OutputDesign', num2str(design),'Set', num2str(arg1),'.mat');
        load(flnm);
        for i_batch = 1:N
            NRE.crude(i_batch) =NRE.crude(i_batch)+ ...
                norm(output.true_gamma-output.crude{i_batch}.mean_gamma)/norm(output.true_gamma)/n_rep;
            [~,~,~,temp] = perfcurve(output.local_connected,output.crude{i_batch}.mean_gamma,1);
            AUC.crude(i_batch) =AUC.crude(i_batch)+ temp/n_rep;
            NRE_mu.crude(i_batch) = NRE_mu.crude(i_batch)+ ...
                norm(output.true_mu -output.crude{i_batch}.mean_mu(output.local_connected)')/norm(output.true_mu)/n_rep;
        end
        for ind_t = 1:length(t_seq)
            %             NRE.Gibbs(i_batch)= NRE.Gibbs(i_batch)+ ...
            %                 norm(output.true_gamma-output.Gibbs{i_batch}.mean_gamma')/norm(output.true_gamma)/n_rep;
            NRE.EM(ind_t)= NRE.EM(ind_t)+ ...
                norm(output.true_gamma- output.EM{ind_t}.mean_gamma)/norm(output.true_gamma)/n_rep;
            
            %             [~,~,~,temp] = perfcurve(output.local_connected,output.Gibbs{i_batch}.mean_gamma,1);
            %             AUC.Gibbs(i_batch) =AUC.Gibbs(i_batch)+ temp/n_rep;
            [~,~,~,temp] = perfcurve(output.local_connected,output.EM{ind_t}.mean_gamma,1);
            AUC.EM(ind_t) =AUC.EM(ind_t) + temp/n_rep;
            
            %             NRE_mu.Gibbs(i_batch) =NRE_mu.Gibbs(i_batch) + ...
            %                 norm(output.true_mu -output.Gibbs{i_batch}.mean_mu(output.local_connected)')/norm(output.true_mu)/n_rep;
            NRE_mu.EM(ind_t) = NRE_mu.EM(ind_t) + ...
                norm(output.true_mu -output.EM{ind_t}.mean_mu(output.local_connected))/norm(output.true_mu)/n_rep;
        end
    end
    summary_stats{design+1}.NRE=NRE;
    summary_stats{design+1}.AUC=AUC;
    summary_stats{design+1}.NRE_mu=NRE_mu;
end

%% Plotting
% figure(1)
% plot(1:N, mean(dt_random_sum,2) ,'col',[1,0,0,1],'Linewidth',4);
% hold on;
% plot(1:N, mean(dt_optimal_sum,2) ,'col',[0,0,1,1],'Linewidth',4);
% ylim([0,8]);
% xlim([0,N]);
% for i = 1:num_sim
%     plot(1:N, dt_random_sum(:,i),'col',[1,0,0,0.1],'Linewidth',1);
%     plot(1:N, dt_optimal_sum(:,i),'col',[0,0,1,0.1],'Linewidth',1);
% end
% 
% xlabel('Number of batches');
% ylabel('Computing time per batch (seconds)');
% hold off;

% NRE of connectivity reconstruction
figure(2)
plot(1:N,summary_stats{1}.NRE.crude,'col',[1,0,0,1],'Linewidth',1);
hold on;
plot(t_seq, summary_stats{1}.NRE.EM,'col',[0,1,0,1],'Linewidth',1);
hold on;
% plot(1:N, summary_stats{1}.NRE.Gibbs,'col',[0,0,1,1],'Linewidth',1);
% hold on;

plot(1:N,summary_stats{2}.NRE.crude,':','col',[1,0,0,1],'Linewidth',2);
hold on;
plot(t_seq, summary_stats{2}.NRE.EM,':','col',[0,1,0,1],'Linewidth',2);
hold on;
% plot(1:N, summary_stats{2}.NRE.Gibbs,':','col',[0,0,1,1],'Linewidth',1);
% hold on;

xlim([0,N+1]);
ylim([0,3]);
xlabel('Number of batches');
ylabel('NRE of connectivity');
%xticks([0 50 100 150 200])
%xticklabels({'0', '50', '100', '150', '200'})
hold off;



figure(3)
plot(1:N,summary_stats{1}.NRE_mu.crude,'col',[1,0,0,1],'Linewidth',1);
hold on;
plot(t_seq, summary_stats{1}.NRE_mu.EM,'col',[0,1,0,1],'Linewidth',1);
hold on;
% plot(1:N, summary_stats{1}.NRE_mu.Gibbs,'col',[0,0,1,1],'Linewidth',1);
% hold on;

plot(1:N,summary_stats{2}.NRE_mu.crude,':','col',[1,0,0,1],'Linewidth',2);
hold on;
plot(t_seq, summary_stats{2}.NRE_mu.EM,':','col',[0,1,0,1],'Linewidth',2);
hold on;
% plot(1:N, summary_stats{2}.NRE_mu.Gibbs,':','col',[0,0,1,1],'Linewidth',1);
% hold on;
xlim([0,N+1]);
ylim([0,5]);
xlabel('Number of batches');
ylabel('NRE of event sizes');
%xticks([0 50 100 150 200])
%xticklabels({'0', '50', '100', '150', '200'})
hold off;


figure(4)
plot(1:N,summary_stats{1}.AUC.crude,'col',[1,0,0,1],'Linewidth',1);
hold on;
plot(t_seq, summary_stats{1}.AUC.EM,'col',[0,1,0,1],'Linewidth',1);
hold on;
% plot(1:N, summary_stats{1}.AUC.Gibbs,'col',[0,0,1,1],'Linewidth',1);
% hold on;

plot(1:N,summary_stats{2}.AUC.crude,':','col',[1,0,0,1],'Linewidth',2);
hold on;
plot(t_seq, summary_stats{2}.AUC.EM,':','col',[0,1,0,1],'Linewidth',2);
hold on;
% plot(1:N, summary_stats{2}.AUC.Gibbs,':','col',[0,0,1,1],'Linewidth',1);
% hold on;

xlim([0,N+1]);
ylim([0,1]);
xlabel('Number of batches');
ylabel('AUC of connectivity');
%xticks([0 50 100 150 200])
%xticklabels({'0', '50', '100', '150', '200'})
hold off;

% saveas(1,strcat(outflnm,'Time','.jpg'));
% 
% saveas(2,strcat(outflnm,'NRE_conn','.jpg'));
% 
% saveas(3,strcat(outflnm,'NRE_mark','.jpg'));
% 
% saveas(4,strcat(outflnm,'AUC_conn','.jpg'));
