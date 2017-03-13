%% Read the estimates from the working model
 flnm=strcat('./Data/Full_crude.mat');
 load(flnm)
mean_gamma_crude = overall_connectivity;
mean_mu_crude = overall_mark;
comptime_crude= t_delta;
%% Read the estimates from the EM algorithm
 flnm=strcat('./Data/Full_EM.mat');
 load(flnm)
mean_gamma_EM = gamma_samples(:,end);
mean_mu_EM = mu_samples(:,end);
mean_sigma_EM= sigma_samples(:,end);
comptime_EM=t_delta;

%% Read the estimates from the Gibbs sampler (with soft assignments)
 flnm=strcat('./Data/Full_minibatch_int.mat');
 load(flnm)
mean_gamma_mini_gibbs =  mean(gamma_samples, 1);
mean_mu_mini_gibbs =  mean(mu_samples, 1);
mean_sigma_mini_gibbs =  mean(sigma_samples, 1);
comptime_mini_gibbs = t_delta;

mean_gamma_mini_gibbs(local_connected)
%% Read the truth of the simulation
flnm=strcat('./Data/truth.mat');
load(flnm);

local_connected =local_neuron_amplitudes>0;
n_connected = sum(local_connected);
n_disconnected = n_cell_local - n_connected;
%% Calculate summary statistics:
% Normalized reconstruction error of connectivity
true_gamma =(local_neuron_amplitudes>0)*(1-evoked_params.failure_prob);
NRE_conn = norm(true_gamma-mean_gamma_crude)/norm(true_gamma);
NRE_conn_mini_gibbs=  norm(true_gamma-mean_gamma_mini_gibbs')/norm(true_gamma);
NRE_conn_EM=  norm(true_gamma-mean_gamma_EM)/norm(true_gamma);

[~,~,~,temp] = perfcurve(local_neuron_amplitudes>0,overall_connectivity ,1);
AUC_conn = temp;


true_mu = local_neuron_amplitudes(local_connected);
NRE_mark = norm(true_mu -mean_mu_crude(local_connected)')/norm(true_mu);
NRE_mark_mini_gibbs = norm(true_mu -mean_mu_mini_gibbs(local_connected)')/norm(true_mu);
NRE_mark_EM = norm(true_mu -mean_mu_EM(local_connected))/norm(true_mu);

%%
% True mean amplitudes v.s. estimates
figure(1)
% With sigma assumed to be known and minibatch
plot(true_mu,mean_mu_mini_gibbs(local_connected)', '.','markers',20,'col',[0,0,1,0.8]);
hold on;
plot(true_mu,mean_mu_EM(local_connected)', '.','markers',20,'col',[0,1,0,0.8]);
hold on;
plot(true_mu,mean_mu_crude(local_connected)', '.','markers',20,'col',[1,0,0,0.8]);
line([0 9],[0 9]); 
ylim([0,max(true_mu)]);
xlim([0,max(true_mu)]);
xlabel('True mean amplitudes');
ylabel('Estimates');
hold off;



figure(2)
% With sigma assumed to be known and minibatch
plot(true_gamma+normrnd(0, 0.04,[n_cell_local 1]),mean_gamma_mini_gibbs', '.','markers',20,'col',[0,0,1,0.8]);
hold on;
plot(true_gamma+normrnd(0, 0.04,[n_cell_local 1]),mean_gamma_EM', '.','markers',20,'col',[0,1,0,0.8]);
hold on;
plot(true_gamma+normrnd(0, 0.04,[n_cell_local 1]),mean_gamma_crude', '.','markers',20,'col',[1,0,0,0.8]);
line([0 1],[0 1]); 
ylim([-0.1,1.1]);
xlim([-0.1,1.1]);
xlabel('True synaptic success rates');
ylabel('Estimates');
hold off;

flnm=strcat('../../Figures/Full_model/');
saveas(1,strcat(flnm,'Gamma','.jpg'));
saveas(2,strcat(flnm,'mu','.jpg'));

%% 
figure(3)
% With sigma assumed to be known and minibatch
plot(true_mu,mean_gamma_mini_gibbs(local_connected)', '.','markers',20,'col',[0,0,1,0.8]);
hold on;
plot(true_mu,mean_gamma_EM(local_connected)', '.','markers',20,'col',[0,1,0,0.8]);
hold on;
plot(true_mu,mean_gamma_crude(local_connected)', '.','markers',20,'col',[1,0,0,0.8]);
ylim([-0.1,1.1]);
xlim([0,max(true_mu)]);
xlabel('True amplitudes');
ylabel('Estimated synaptic success rate');
hold off;
saveas(3,strcat(flnm,'GammavsMu','.jpg'));

    %%
     NRE_conn
     NRE_conn_mini_gibbs
     NRE_conn_EM
      
     NRE_mark
     NRE_mark_mini_gibbs
     NRE_mark_EM

     %% Variance:
     mean(mean_sigma_EM(local_connected))
     std(mean_sigma_EM(local_connected))
     
    mean(mean_sigma_mini_gibbs(local_connected))
    std(mean_sigma_mini_gibbs(local_connected))
    
   

%%
comptime_EM
comptime_crude
comptime_mini_gibbs

%% Fake cells v.s. true cells 
Z(local_connected,1:2)

dis_measure = zeros(n_cell_local,1);

for i = 1:n_cell_local
    total_dist = ones(sum(local_connected),1)*Z(i,1:2)-Z(local_connected,1:2);
    sq_disc = total_dist.^2; 
    dis_measure(i) = min(sqrt(sum(sq_disc,2)));
end

%% Estimated synaptic sucess rate vs distance to true connected cells
figure(4)
% With sigma assumed to be known and minibatch
plot(dis_measure,mean_gamma_mini_gibbs', '.','markers',20,'col',[0,0,1,0.8]);
hold on;
plot(dis_measure,mean_gamma_EM', '.','markers',20,'col',[0,1,0,0.8]);
hold on;
plot(dis_measure,mean_gamma_crude', '.','markers',20,'col',[1,0,0,0.8]);
ylim([-0.1,1.1]);
xlim([-5,max(dis_measure)]);
xlabel('Distance to the closest connected cell');
ylabel('Estimated synaptic success rate');
hold off;
saveas(4,strcat(flnm,'GammavsDist','.jpg'));
