%% Evaluate the performance:
% 1, normalized reconstruction error of connectivity
% 2, AUC for connectivity
% 3, normalized amplitudes reconstruction
% Z_full = [local_neuron_amplitudes>0 Z];
% density_true = ksdensity2d(Z_full,x_dense,y_dense,bandwidth);
% Z_amp = [local_neuron_amplitudes Z];
% amp_true = ksdensity2d(Z_amp,x_dense,y_dense,bandwidth);


NRE_conn = zeros(N/10,1);
NRE_amp = zeros(N/10,1);
AUC_conn = zeros(N/10,1);

Spatial_conn = zeros(N/10,1);
Spatial_amp = zeros(N/10,1);

for outer_i = 1:(N/10)
    overall_connectivity = zeros(K_z,1);
    overall_connectivity_w = zeros(K_z,1);
    
    overall_amplitudes = zeros(K_z,1);
    normalized_constants = zeros(K_z,1);
    for j = 1:num_threshold
        w_estimate = output_eva{outer_i}(j).mu(2:end).*output_eva{outer_i}(j).alpha(2:end);
        %w_estimate = output_eva{outer_i}(j).alpha(2:end);
        
        overall_connectivity = max(overall_connectivity, ...
     output_eva{outer_i}(j).alpha(2:end));
        overall_connectivity_w = overall_connectivity_w+  w_estimate;
        overall_amplitudes = overall_amplitudes + ...
            (amplitude_threshold(j)+amplitude_threshold(j+1)) *w_estimate/2;
    end
    
    % Normalized reconstruction error of connectivity
    NRE_conn(outer_i) = norm((local_neuron_amplitudes>0) -overall_connectivity_w )/norm( 1*(local_neuron_amplitudes(:,1)>0));
    NRE_amp(outer_i) = norm(local_neuron_amplitudes -overall_amplitudes)/norm(local_neuron_amplitudes(:,1));
    
    [~,~,~,temp] = perfcurve(local_neuron_amplitudes>0,overall_connectivity_w ,1);
    AUC_conn(outer_i) = temp;
    
%     Z_est = [overall_connectivity_w Z];
%     density_est = ksdensity2d(Z_est,x_dense,y_dense,bandwidth);
%        
%     Z_ampest = [overall_amplitudes Z];
%     amp_est = ksdensity2d(Z_ampest,x_dense,y_dense,bandwidth);
%     
%     Spatial_amp(outer_i) = norm(amp_true-amp_est,2);
%     Spatial_conn(outer_i) = norm(density_true-density_est,2);
    
end
