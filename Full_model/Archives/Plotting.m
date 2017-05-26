temp_struct = struct;
temp_struct.crude=0;
temp_struct.Gibbs=0;
temp_struct.EM=0;

 NRE=temp_struct;
 AUC=temp_struct;
 NRE_mu= temp_struct;
 
 n_rep=50;
for arg1 = 1:n_rep
flnm=strcat('./Results/Output', num2str(arg1),'.mat');

load(flnm);

NRE.crude =NRE.crude+ norm(output.true_gamma-output.crude.mean_gamma)/norm(output.true_gamma)/n_rep;
NRE.Gibbs= NRE.Gibbs+ norm(output.true_gamma-output.Gibbs.mean_gamma')/norm(output.true_gamma)/n_rep;
NRE.EM= NRE.EM+ norm(output.true_gamma- output.EM.mean_gamma)/norm(output.true_gamma)/n_rep;

[~,~,~,temp] = perfcurve(output.local_connected,output.crude.mean_gamma,1);
AUC.crude =AUC.crude+ temp/n_rep;
[~,~,~,temp] = perfcurve(output.local_connected,output.Gibbs.mean_gamma,1);
AUC.Gibbs =AUC.Gibbs+ temp/n_rep;
[~,~,~,temp] = perfcurve(output.local_connected,output.EM.mean_gamma,1);
AUC.EM =AUC.EM + temp/n_rep;



NRE_mu.crude = NRE_mu.crude+norm(output.true_mu -output.crude.mean_mu(output.local_connected)')/norm(output.true_mu)/n_rep;
NRE_mu.Gibbs =NRE_mu.Gibbs + norm(output.true_mu -output.Gibbs.mean_mu(output.local_connected)')/norm(output.true_mu)/n_rep;
NRE_mu.EM = NRE_mu.EM +norm(output.true_mu -output.EM.mean_mu(output.local_connected))/norm(output.true_mu)/n_rep;
end
%%
NRE
AUC
NRE_mu


%% Collapse the estimates for plotting
temp_struct = struct;
temp_struct.mu=[];
temp_struct.gamma=[];
temp_struct.sigma=[];

true_value = temp_struct;
true_value.dis_measure=[]; 
EM=temp_struct;
 Gibbs=temp_struct;
 crude= temp_struct;
 
 
 %n_rep=1;
 arg1 = 3;
%for arg1 = 1:n_rep
flnm=strcat('./Results/Output', num2str(arg1),'.mat');

load(flnm);

true_value.gamma = [true_value.gamma; output.true_gamma];

EM.gamma = [EM.gamma; output.EM.mean_gamma];
crude.gamma = [crude.gamma; output.crude.mean_gamma];
Gibbs.gamma = [Gibbs.gamma; output.Gibbs.mean_gamma'];


true_value.mu = [true_value.mu; output.true_mu];
EM.mu = [EM.mu; output.EM.mean_mu];
crude.mu = [crude.mu; output.crude.mean_mu];
Gibbs.mu = [Gibbs.mu; output.Gibbs.mean_mu'];

EM.sigma = [EM.sigma; output.EM.mean_sigma];
Gibbs.sigma = [Gibbs.sigma; output.Gibbs.mean_sigma'];

true_value.dis_measure = [true_value.dis_measure; output.dis_measure];

%end

%%

% True mean amplitudes v.s. estimates
figure(1)
% With sigma assumed to be known and minibatch
hg=plot(true_value.mu,Gibbs.mu(true_value.gamma>0), '.','markers',20,'col',[0,0,1,0.4]);
hold on;
he=plot(true_value.mu,EM.mu(true_value.gamma>0), '.','markers',20,'col',[0,1,0,0.4]);
hold on;
hw=plot(true_value.mu,crude.mu(true_value.gamma>0), '.','markers',20,'col',[1,0,0,0.4]);
line([0 9],[0 9]); 
ylim([0,max(true_value.mu)]);
xlim([0,max(true_value.mu)]);
xlabel('True mean amplitudes');
ylabel('Estimates');
legend([hw,hg,he],'Working', 'Gibbs', 'EM','Location','southeast');

hold off;

%{[1,0,0,0.4],[0,0,1,0.4],[0,1,0,0.4] }

figure(2)
hg=plot(true_value.gamma+normrnd(0, 0.04,[size(true_value.gamma,1), 1]),Gibbs.gamma, '.','markers',20,'col',[0,0,1,0.8]);
hold on;
he=plot(true_value.gamma+normrnd(0, 0.04,[size(true_value.gamma,1), 1]),EM.gamma, '.','markers',20,'col',[0,1,0,0.8]);
hold on;
hw=plot(true_value.gamma+normrnd(0, 0.04,[size(true_value.gamma,1), 1]),crude.gamma, '.','markers',20,'col',[1,0,0,0.8]);
line([0 1],[0 1]); 
ylim([-0.1,1.1]);
xlim([-0.1,1.1]);
xlabel('True synaptic success rates');
ylabel('Estimates');
legend([hw,hg,he],'Working', 'Gibbs', 'EM','Location','southeast');

hold off;

flnm=strcat('../../Figures/Full_model/');
saveas(1,strcat(flnm,'Gamma','.jpg'));
saveas(2,strcat(flnm,'mu','.jpg'));

%% Estimated synaptic sucess rate vs distance to true connected cells
figure(4)
% With sigma assumed to be known and minibatch
hg=plot(true_value.dis_measure,Gibbs.gamma, '.','markers',20,'col',[0,0,1,0.4]);
hold on;
he=plot(true_value.dis_measure,EM.gamma, '.','markers',20,'col',[0,1,0,0.4]);
hold on;
hw= plot(true_value.dis_measure,crude.gamma, '.','markers',20,'col',[1,0,0,0.4]);
ylim([-0.1,1.1]);
xlim([-5,max(true_value.dis_measure)]);
xlabel('Distance to the closest connected cell');
ylabel('Estimated synaptic success rate');

legend([hw,hg,he],'Working', 'Gibbs', 'EM','Location','northeast');

hold off;
saveas(4,strcat(flnm,'GammavsDist','.jpg'));

%%
mean(Gibbs.sigma(true_value.gamma>0))
std(Gibbs.sigma(true_value.gamma>0))

mean(EM.sigma(true_value.gamma>0))
std(EM.sigma(true_value.gamma>0))