for arg1 = 1:10
    
   
% Read the estimates from the working model
 flnm=strcat('./Data/Full_crude', num2str(arg1),'.mat');
 load(flnm)
crude = struct;
crude.mean_gamma = overall_connectivity;
crude.mean_mu = overall_mark;
crude.comptime= t_delta;
% Read the estimates from the EM algorithm
 flnm=strcat('./Data/Full_EM', num2str(arg1),'.mat');
 load(flnm)
EM.mean_gamma = gamma_samples(:,end);
EM.mean_mu = mu_samples(:,end);
EM.mean_sigma= sigma_samples(:,end);
EM.comptime=t_delta;

% Read the estimates from the Gibbs sampler (with soft assignments)
 flnm=strcat('./Data/Full_minibatch_int', num2str(arg1),'.mat');
 load(flnm)
Gibbs.mean_gamma =  mean(gamma_samples, 1);
Gibbs.mean_mu =  mean(mu_samples, 1);
Gibbs.mean_sigma =  mean(sigma_samples, 1);
Gibbs.comptime = t_delta;


% Read the truth of the simulation
flnm=strcat('./Data/truth', num2str(arg1),'.mat');
load(flnm);

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
true_gamma =(local_amplitudes>0)*(1-evoked_params.failure_prob);
true_mu = local_amplitudes(local_connected);


output.dis_measure = dis_measure;
output.local_amplitudes = local_amplitudes;

output.local_connected = local_connected;
output.true_gamma = true_gamma;
output.true_mu = true_mu;
output.EM= EM;
output.Gibbs = Gibbs;
output.crude=crude;

flnm=strcat('./Results/Output', num2str(arg1),'.mat');

save(flnm,'output');

end
