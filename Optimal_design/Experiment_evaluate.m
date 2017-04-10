%% Evaluate the performance:
% need to turn this into a function too
% 1, normalized reconstruction error of connectivity
% 2, AUC for connectivity
% 3, normalized mark reconstruction (latency or amplitude)

NRE_conn = zeros(N/10,1);
NRE_mark = zeros(N/10,1);
AUC_conn = zeros(N/10,1);

for outer_i = 1:(N/10)
    overall_connectivity = zeros(K_z,1);
    overall_connectivity_w = zeros(K_z,1);
    
    overall_mark = zeros(K_z,1);
    normalized_constants = zeros(K_z,1);
    for j = 1:num_threshold
        w_estimate = output_eva{outer_i}(j).mu(2:end).*output_eva{outer_i}(j).alpha(2:end);
        overall_connectivity = max(overall_connectivity, ...
     output_eva{outer_i}(j).alpha(2:end));
        overall_connectivity_w = overall_connectivity_w+  w_estimate;
        overall_mark = overall_mark + ...
            (threshold(j)+threshold(j+1)) *w_estimate/2;
    end
    
    % Normalized reconstruction error of connectivity
    NRE_conn(outer_i) = norm((local_neuron_amplitudes>0) -overall_connectivity_w )/norm( 1*(local_neuron_amplitudes(:,1)>0));
    if mark == 0
        NRE_mark(outer_i) = norm(local_neuron_amplitudes -overall_mark)/norm(local_neuron_amplitudes(:,1));
    elseif mark == 1
        NRE_mark(outer_i) = norm(local_neuron_latencies/data_params.dt -overall_mark)/norm(local_neuron_amplitudes(:,1));
    end
    
    [~,~,~,temp] = perfcurve(local_neuron_amplitudes>0,overall_connectivity_w,1);
    AUC_conn(outer_i) = temp;
 
    
end
