
for i_setting = [1 2 4:9]
N=200;

t_seq = [1 2 3 4 5 5*(1: N/5)];
   
temp_struct = struct;
temp_struct.crude=zeros(N,1);
temp_struct.EM=zeros(length(t_seq),1);

temp_struct.working=zeros(length(t_seq),1);
summary_stats = cell(2);
i_rep =  [1:50];
n_rep=length(i_rep);

for design = 0:1
    NRE=temp_struct;
    AUC=temp_struct;
    NRE_mu=temp_struct;
    Comp_time = zeros(1,N);
    for arg1 = i_rep
        flnm=strcat('./Results/Apr19/OutputDesign', num2str(design),'Setting',num2str(i_setting),...
                'Rep', num2str(arg1),'.mat');
        load(flnm);
        for i_batch = 1:N
            NRE.crude(i_batch) =NRE.crude(i_batch)+ ...
                norm(output.true_gamma-output.crude{i_batch}.mean_gamma)/norm(output.true_gamma)/n_rep;
            [~,~,~,temp] = perfcurve(output.local_connected,output.crude{i_batch}.mean_gamma,1);
            AUC.crude(i_batch) =AUC.crude(i_batch)+ temp/n_rep;
            NRE_mu.crude(i_batch) = NRE_mu.crude(i_batch)+ ...
                norm(output.true_mu -output.crude{i_batch}.mean_mu(output.local_connected)')/norm(output.true_mu)/n_rep;
           
        end
       Comp_time = Comp_time+  output.crude{i_batch}.comptime(1:N)/n_rep;
        
%         for ind_t = 1:length(t_seq)
%             NRE.working(ind_t)= NRE.working(ind_t)+ ...
%                 norm(output.true_gamma- output.working{ind_t}.mean_gamma)/norm(output.true_gamma)/n_rep;
%             
%             [~,~,~,temp] = perfcurve(output.local_connected,output.working{ind_t}.mean_gamma,1);
%             AUC.working(ind_t) =AUC.working(ind_t) + temp/n_rep;
%             
%             NRE_mu.working(ind_t) = NRE_mu.working(ind_t) + ...
%                 norm(output.true_mu -output.working{ind_t}.mean_mu(output.local_connected))/norm(output.true_mu)/n_rep;
%         end
        
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
    summary_stats{design+1}.comp=Comp_time;
    
end

% Plotting
figure(1)
rw=plot(1:N, summary_stats{1}.comp,'col',[1,0,0,1],'Linewidth',1);
hold on;
ow=plot(1:N, summary_stats{2}.comp,':','col',[1,0,0,1],'Linewidth',2);
hold on;
%legend([rw,ow],'Working - Random', 'Working - Optimal', 'Location','northeast');
ylim([0,3.2]);
xlim([0,N]);

xlabel('Number of batches');
ylabel('Computing time per batch (seconds)');
hold off;

% NRE of connectivity reconstruction
figure(2)
rw=plot(1:N,summary_stats{1}.NRE.crude(1:N),'col',[1,0,0,1],'Linewidth',1);
hold on;
re=plot(t_seq, summary_stats{1}.NRE.EM,'col',[0,1,0,1],'Linewidth',1);
hold on;
ow=plot(1:N,summary_stats{2}.NRE.crude(1:N),':','col',[1,0,0,1],'Linewidth',2);
hold on;
oe=plot(t_seq, summary_stats{2}.NRE.EM,':','col',[0,1,0,1],'Linewidth',2);
hold on;


legend([rw,re,ow,oe],'Working - Random', 'EM - Random', 'Working - Optimal', 'EM - Optimal', 'Location','northeast');

xlim([0,N+1]);
ylim([0,2]);
xlabel('Number of batches');
ylabel('NRE of synaptic success rates');
%xticks([0 50 100 150 200])
%xticklabels({'0', '50', '100', '150', '200'})
hold off;



figure(3)
rw= plot(1:N,summary_stats{1}.NRE_mu.crude(1:N),'col',[1,0,0,1],'Linewidth',1);
hold on;
re= plot(t_seq, summary_stats{1}.NRE_mu.EM,'col',[0,1,0,1],'Linewidth',1);
hold on;

ow=plot(1:N,summary_stats{2}.NRE_mu.crude(1:N),':','col',[1,0,0,1],'Linewidth',2);
hold on;
oe=plot(t_seq, summary_stats{2}.NRE_mu.EM,':','col',[0,1,0,1],'Linewidth',2);
hold on;
legend([rw,re,ow,oe],'Working - Random', 'EM - Random', 'Working - Optimal', 'EM - Optimal', 'Location','northeast');
xlim([0,N+1]);
ylim([0,4]);
xlabel('Number of batches');
ylabel('NRE of mean event sizes');
hold off;


figure(4)
rw=plot(1:N,summary_stats{1}.AUC.crude(1:N),'col',[1,0,0,1],'Linewidth',1);
hold on;
re=plot(t_seq, summary_stats{1}.AUC.EM,'col',[0,1,0,1],'Linewidth',1);
hold on;

ow=plot(1:N,summary_stats{2}.AUC.crude(1:N),':','col',[1,0,0,1],'Linewidth',2);
hold on;
oe=plot(t_seq, summary_stats{2}.AUC.EM,':','col',[0,1,0,1],'Linewidth',2);
hold on;

legend([rw,re,ow,oe],'Working - Random', 'EM - Random', 'Working - Optimal', 'EM - Optimal', 'Location','southeast');
xlim([0,N+1]);
ylim([0.5,1]);
xlabel('Number of batches');
ylabel('AUC');
hold off;

outflnm =strcat('./Figures/Apr19/PlotSetting',num2str(i_setting)); 

saveas(1,strcat(outflnm,'time','.jpg'));
saveas(2,strcat(outflnm,'NRE_conn','.jpg'));
saveas(3,strcat(outflnm,'NRE_mark','.jpg'));
saveas(4,strcat(outflnm,'AUC_conn','.jpg'));
end
