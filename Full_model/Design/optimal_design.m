function [locations_trials,counts_freq] = optimal_design(...
    i_batch, num_sources,num_peaks,num_trials_first,num_trials_batch, output, Y_g, ...
    random_prop, counts_freq, pi_dense_local, inner_normalized_products,grid_index, freq_pen, num_samples)
fprintf('Batch %d\n', i_batch);
% Designing trials
if i_batch  == 1
    num_this_batch = num_trials_first;
    counts_freq =  zeros(size(pi_dense_local,2),1);
    ind_seq = zeros(0,1);
    locations_trials = zeros(num_this_batch ,num_sources);
    num_sequence  = ceil(num_this_batch*num_sources/length(grid_index));
    temp_index = 1:length(grid_index);
    for i_seq = 1:num_sequence 
        ind_seq = [ind_seq randsample(temp_index, length(temp_index),false)]; 
    end
    for l = 1:num_this_batch
        locations_next = ind_seq(  (1: (num_sources))+ (l-1)*(num_sources) );
        locations_trials(l,:)=locations_next;
        counts_freq(locations_next) = counts_freq(locations_next)+1;
    end
else 
    num_this_batch = num_trials_batch;
    locations_trials = zeros(num_this_batch ,num_sources);
    
    % Optiml design
    if random_prop < 1
        %--------------------------------------------------%
        % Part I: evaluate the entropies
        n_cell_local = length(output(1).alpha)-1;
        num_threshold = length(output);
        
        delta_H = zeros(n_cell_local+1,1);
        
        for j = 1:(num_threshold)
            H_current = per_neuron_joint_entropy(output(j).alpha, ...
                output(j).mu, output(j).s_sq);
            Y_n = Y_g(:,j);
            if  sum(Y_n)==0 % prevent zero standard deviation
                hyperparam_sigma_n = 1;
            else
                hyperparam_sigma_n = std(Y_n);
            end
            H_expected = approximate_expected_joint_entropy_single_neuron_poisson(...
                output(j).alpha, output(j).mu, ...
                output(j).s_sq, hyperparam_sigma_n, num_samples);
            delta_H = delta_H + H_current-H_expected;
        end
        entropy_locations = pi_dense_local'*delta_H(2:end);
        % Encourage stimulating rare locations (e.g. nearby cells)
        weights_freq= freq_pen -counts_freq/max(counts_freq);
        % use 1.1 so that the penalty will not be 0.
        entropy_locations = weights_freq.*entropy_locations;
        
        %-----------------------------------------------------------%
        % Part II: find peaks and random locations 
        idx = [];
        for i = 1:num_peaks
            [m_id temp_idx] = max(entropy_locations);
            idx = [idx temp_idx];
            entropy_locations = entropy_locations - inner_normalized_products(:,temp_idx)*m_id;
        end
    end
    
    %----------------------------------------------------------%
    % Part III: Determine the locations 
    % Find multiple treaments:
    % Permute the index vector
    num_rep = ceil(num_sources*num_this_batch/num_peaks);
    idx_seq = zeros(num_rep*num_peaks,1);
    for l = 1:num_rep
        idx_seq((l-1)*num_peaks+(1:num_peaks)) = randsample(idx, num_peaks);
    end
    % Replace some locations with random locations 
    num_random_peaks = floor(random_prop*length(idx_seq));
    idx_seq(randsample(1:length(idx_seq), num_random_peaks)) = randsample(grid_index, num_random_peaks);
    for l = 1:num_this_batch
        locations_next = idx_seq((l-1)*num_sources + (1:num_sources));
        locations_trials(l,:)=locations_next;
        counts_freq(locations_next) = counts_freq(locations_next)+1;
    end
end