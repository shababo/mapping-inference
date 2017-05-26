function [locations_trials,powers_trials,counts_freq] = optimal_design_v2(...
    i_batch, num_sources,num_peaks,num_trials_first,num_trials_batch, output, Y_g, ...
    num_power_level,...
    random_prop, counts_freq, pi_dense_local, inner_normalized_products,grid_index, freq_pen, num_samples)
fprintf('Batch %d\n', i_batch);
% Designing trials
if i_batch  == 1
    num_this_batch = num_trials_first;
    
    ind_seq = zeros(0,1);
    locations_trials = zeros(num_this_batch ,num_sources);
    
    powers_trials= zeros(num_this_batch ,num_sources);
    num_sequence  = ceil(num_this_batch*num_sources/length(grid_index));
    temp_index = 1:length(grid_index);
    for i_seq = 1:num_sequence 
        ind_seq = [ind_seq randsample(temp_index, length(temp_index),false)]; 
    end
    for l = 1:num_this_batch
        locations_next = ind_seq(  (1: (num_sources))+ (l-1)*(num_sources) );
        locations_trials(l,:)=locations_next;
        powers_trials(l,:)= randsample(1:num_power_level, num_sources,true);
        for m= 1:num_sources
                idx_vec = powers_trials(l,m)+ (locations_next(m)-1)*num_power_level;
                counts_freq(idx_vec)=counts_freq(idx_vec)+1;
        end
    end
else 
    num_this_batch = num_trials_batch;
    locations_trials = zeros(num_this_batch ,num_sources);
    powers_trials = zeros(num_this_batch ,num_sources);
    n_grid = size(pi_dense_local,2);
    % Optiml design
    if random_prop < 1
        %--------------------------------------------------%
        % Part I: evaluate the entropies
        n_cell_local = (length(output(1).alpha)-1)/num_power_level;
        num_threshold = length(output);
        
        delta_H = zeros(n_cell_local*num_power_level+1,1);
        
        % Need to modify this part to allow for different powers:
        
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
        
        % Calculate the entropy at each location for each power level:
        entropy_locations=zeros(n_grid*num_power_level,1);
        weights_freq=zeros(n_grid,1);
        for i_pl = 1:num_power_level
            ind_cell_vec = 1+i_pl+ (0:n_cell_local-1)*num_power_level;
            ind_grid_vec = i_pl+ (0:n_grid-1)*num_power_level;
            
            entropy_locations(ind_grid_vec) = pi_dense_local'*delta_H(ind_cell_vec);
            % Encourage stimulating rare locations (e.g. nearby cells)
            weights_freq= freq_pen -counts_freq( ind_grid_vec )/max(counts_freq);
            entropy_locations(ind_grid_vec) = weights_freq.*entropy_locations(ind_grid_vec);
        end
        
        % use 1.1 so that the penalty will not be 0.
        
        %-----------------------------------------------------------%
        % Part II: find peaks and random locations 
        idx = [];
        idx_mat = zeros(n_cell_local,num_power_level);
        
        for i = 1:num_peaks
            [m_id temp_idx] = max(entropy_locations);
            idx = [idx temp_idx];
            i_pl = 1+mod(temp_idx-1, num_power_level); 
            idx_loca = (temp_idx - i_pl)/num_power_level + 1;
            % Find out which power level the peak belongs to 
             ind_grid_vec = i_pl+ (0:n_grid-1)*num_power_level;
            entropy_locations(ind_grid_vec) = entropy_locations(ind_grid_vec) - inner_normalized_products(:,idx_loca)*m_id;
            idx_mat(idx_loca,i_pl)=1;
        end
        
        idx_locations =find( sum(idx_mat,2));
        
        
    %----------------------------------------------------------%
    % Part III: Determine the locations 
    % Find multiple treaments:
    % Permute the index vector
    num_rep = ceil(num_sources*num_this_batch/length(idx_locations));
    idx_seq = zeros(num_rep*length(idx_locations),1);
    for l = 1:num_rep
        idx_seq((l-1)*length(idx_locations) +(1:length(idx_locations))) = randsample(idx_locations, length(idx_locations));
    end
    % Replace some locations with random locations 
    num_random_peaks = floor(random_prop*length(idx_seq));
    idx_not_selected =find( sum(idx_mat,2)==0);
    idx_replacement = randsample(idx_not_selected, num_random_peaks);
    idx_seq(randsample(1:length(idx_seq), num_random_peaks)) = idx_replacement;
    idx_mat(unique(idx_replacement),:)=1; % allow all powers in the random locations 
    for l = 1:num_this_batch
        locations_next = idx_seq((l-1)*num_sources + (1:num_sources));
        locations_trials(l,:)=locations_next;
        for m= 1:num_sources
            idx_power_level = find(idx_mat(locations_next(m),:));
            powers_trials(l,m) = randsample(idx_power_level, 1);
            idx_vec = powers_trials(l,m)+ (locations_next(m)-1)*num_power_level;
            counts_freq(idx_vec)=counts_freq(idx_vec)+1;
        end
    end
    
   
    else % random design 
         ind_seq = zeros(0,1);
   
         num_sequence  = ceil(num_this_batch*num_sources/length(grid_index));
    temp_index = 1:length(grid_index);
    for i_seq = 1:num_sequence 
        ind_seq = [ind_seq randsample(temp_index, length(temp_index),false)]; 
    end
        for l = 1:num_this_batch
            locations_next = ind_seq(  (1: (num_sources))+ (l-1)*(num_sources) );
            locations_trials(l,:)=locations_next;
            powers_trials(l,:)= randsample(1:num_power_level, num_sources,true);
            for m= 1:num_sources
                idx_vec = powers_trials(l,m)+ (locations_next(m)-1)*num_power_level;
                counts_freq(idx_vec)=counts_freq(idx_vec)+1;
            end
        
        end
    end

end


