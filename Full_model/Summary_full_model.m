
 flnm=strcat('./Data/Full_minibatch_int.mat');
 load(flnm)
mean_gamma_known_mini_int = mean(gamma_samples,1);
mean_amp_known_mini_int = mean(mu_samples,1);
mean_sigma_known_mini_int= mean(sigma_samples,1);
comptime_known_mini_int=t_delta;

local_connected =local_neuron_amplitudes>0;
n_connected = sum(local_connected);
n_disconnected = n_cell_local - n_connected;

 flnm=strcat('./Data/Full_minibatch_nomark.mat');
 load(flnm)
mean_gamma_known_mini_nomark = mean(gamma_samples,1);
mean_amp_known_mini_nomark = mean(mu_samples,1);
mean_sigma_known_mini_nomark= mean(sigma_samples,1);
comptime_known_mini_nomark=t_delta;

figure(1)


% With sigma assumed to be known and minibatch
plot(1:n_connected,mean_gamma_known_mini_int(find(local_connected)),'.','markers',20,'col',[0,0,1,0.8]);
hold on;
plot(n_connected+(1:n_disconnected),mean_gamma_known_mini_int(find(local_connected==0)),'.','markers',10,'col',[0,0,1,0.8]);
hold on;


% With sigma assumed to be known and minibatch
plot(1:n_connected,mean_gamma_known_mini_nomark(find(local_connected)),'.','markers',20,'col',[0,1,0,0.8]);
hold on;
plot(n_connected+(1:n_disconnected),mean_gamma_known_mini_nomark(find(local_connected==0)),'.','markers',10,'col',[0,1,0,0.8]);
hold on;

% The crude estimates 
plot(1:n_connected,overall_connectivity(find(local_connected)),'o','markers',10,'col',[1,0,0,0.8]);
hold on;
plot(n_connected+(1:n_disconnected),overall_connectivity(find(local_connected==0)),'o','markers',10,'col',[1,0,0,0.8]);
hold on;
line([0 0],[1 0]); 
line([n_connected+0.5 n_connected+0.5],[8 0]); 
ylim([0,1]);
xlim([0,n_cell_local]);
xlabel('Cells');
ylabel('Posterior mean of Gamma');
hold off;

%% Calculate summary statistics:
  % Normalized reconstruction error of connectivity
    NRE_conn = norm((local_neuron_amplitudes>0)*(1-evoked_params.failure_prob)-overall_connectivity)/...
        norm( (1-evoked_params.failure_prob)*(local_neuron_amplitudes(:,1)>0));
        NRE_mark= norm(local_neuron_amplitudes(local_neuron_amplitudes>0) -overall_mark(local_neuron_amplitudes>0))/...
            norm(local_neuron_amplitudes(local_neuron_amplitudes>0));
    
    [~,~,~,temp] = perfcurve(local_neuron_amplitudes>0,overall_connectivity ,1);
    AUC_conn = temp;
 

 NRE_conn_nomark = norm((local_neuron_amplitudes>0)*(1-evoked_params.failure_prob) -mean_gamma_known_mini_nomark')/...
     norm( (1-evoked_params.failure_prob)*(local_neuron_amplitudes(:,1)>0));
 NRE_mark_nomark= norm(local_neuron_amplitudes(local_neuron_amplitudes>0) -mean_amp_known_mini_nomark(local_neuron_amplitudes>0)' )/...
     norm(local_neuron_amplitudes(local_neuron_amplitudes>0));
    
    [~,~,~,temp] = perfcurve(local_neuron_amplitudes>0,mean_gamma_known_mini_nomark',1);
    AUC_conn_nomark = temp;
 
    
    

 NRE_conn_int = norm((local_neuron_amplitudes>0)*(1-evoked_params.failure_prob) -mean_gamma_known_mini_int')/...
     norm( (1-evoked_params.failure_prob)*(local_neuron_amplitudes(:,1)>0));
 NRE_mark_int= norm(local_neuron_amplitudes(local_neuron_amplitudes>0) -mean_amp_known_mini_int(local_neuron_amplitudes>0)' )/...
     norm(local_neuron_amplitudes(local_neuron_amplitudes>0));
    
    [~,~,~,temp] = perfcurve(local_neuron_amplitudes>0,mean_gamma_known_mini_int',1);
    AUC_conn_int = temp;
 
    %%
     NRE_conn
     NRE_conn_int
     NRE_conn_nomark
      
     NRE_mark
     NRE_mark_int
     NRE_mark_nomark
     
     AUC_conn
     AUC_conn_int
     AUC_conn_nomark
      
