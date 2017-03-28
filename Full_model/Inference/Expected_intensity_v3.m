%%Estimate the intensity for batch data 
%% Parameters:
%

n_trial_temp = size(stimuli_temp,1);
stimuli_bins = cell(size(stimuli_temp,2),1);
M_grid_intensity = cell(size(stimuli_temp,2),n_stimuli_grid);
% marginal expectation
M_intensity=cell(size(stimuli_temp));
n_cell_local = size(stimuli_temp,2);

n_grid_voltage = 200;
n_grid_time = length(t_vect);
t_grid=t_vect;
t_factor = 1;
% Prepare data:
t_grid_upp = t_grid+dt/2;
t_grid_low = t_grid-dt/2;
for i_cell = 1:n_cell_local
    V_th = V_thresholds(i_cell);
    V_reset = V_resets(i_cell);
    E_L=0;
    V_max = V_th+5;
    V_min= V_reset-5;
    v_grid =  (V_max -V_min)*(0: (n_grid_voltage-1) )/(n_grid_voltage-1) +V_min;
    
    [~, index_reset] = min(abs(v_grid-V_reset));
    [~, index_rest] = min(abs(v_grid-E_L));
    
    % Obtain stimuli information on this cell
    if sum(stimuli_temp(:,i_cell)> k_minimum ) > 5
        
        stim_seq =stimuli_temp(stimuli_temp(:,i_cell)> k_minimum,i_cell);
        gap_stimuli=range(stim_seq)/n_stimuli_grid;
        stimuli_bins{i_cell} = min(stim_seq):gap_stimuli:max(stim_seq);
        mid_points = (stimuli_bins{i_cell}(2:end)+stimuli_bins{i_cell}(1:(n_stimuli_grid)))/2;
        
        for i_stimuli = 1:n_stimuli_grid
            k_temp = mid_points(i_stimuli);
            pL_given_V = zeros([2 n_grid_voltage]);
            pL_given_V(2,:) = min(1,t_factor*max(v_grid - V_th,0));
            pL_given_V(1,:) = 1-pL_given_V(2,:);
            pVL_given_I = zeros([n_grid_time n_grid_voltage 2]);
            
            pVL_given_I(1,index_rest,1)=1;
            
            for i_t = 2:n_grid_time
                pVnext_given_V_L = zeros([n_grid_voltage n_grid_voltage 2]);
                pVnext_given_V_L(index_reset,:,2) = 1;
                % We need to calculte this for each time point since the current changes
                % each time
                for i_v = 1:n_grid_voltage
                    v_noise = v_grid-v_grid(i_v)-((E_L-v_grid(i_v))*g + I_stimuli(i_t)*k_temp)*dt;
                    relevant_index = (v_noise > (stoc_mu-3*stoc_sigma)) & (v_noise < (stoc_mu+3*stoc_sigma));
                    pVnext_given_V_L(relevant_index ,i_v,1) = normpdf(v_noise(relevant_index),stoc_mu,stoc_sigma)*(v_grid(2)-v_grid(1));
                end
                %pVnext_given_VL(v,vp,ip)
                %pL_given_V(i,v)
                %pVL_given_I(i_t-1,vp,ip)
                for i_v = 1:n_grid_voltage
                    temp_II_and_III = pVnext_given_V_L(i_v,:,:).*pVL_given_I(i_t-1,:,:);
                    for i = 1:2
                        pVL_given_I(i_t, i_v, i)=pL_given_V(i,i_v)*sum(sum(temp_II_and_III));
                    end
                end
            end
            
            M_grid_intensity{i_cell, i_stimuli} = sum(pVL_given_I(:,:,2),2);
            
        end
        for i_trial = 1:n_trial_temp
            
            % Initialize the storage
            % Simulated traces for marginal intensity
            k = stimuli_temp(i_trial, i_cell);
            if k > k_minimum
                [~, ind_seq] = min(abs(k-mid_points));
                M_intensity{i_trial,i_cell} = M_grid_intensity{i_cell, ind_seq};
            else
                M_intensity{i_trial,i_cell} = zeros([length(t_vect) 1]);
            end
        end
    else % stimulated in fewer than 5 trials
        
        for i_trial = 1:n_trial
            M_intensity{i_trial,i_cell} = zeros([n_grid_time 1]);
            k = stimuli_temp(i_trial, i_cell);
            if k > k_minimum
                pL_given_V = zeros([2 n_grid_voltage]);
                pL_given_V(2,:) = min(1,t_factor*max(v_grid - V_thresholds(i_cell),0));
                pL_given_V(1,:) = 1-pL_given_V(2,:);
                pVL_given_I = zeros([n_grid_time n_grid_voltage 2]);
                pVL_given_I(1,index_rest,1)=1;
                
                %  v_grid-v_grid(i_v)-((E_L-v_grid(i_v))*g + I_stimuli(i_t)*k)*dt
                
                for i_t = 2:n_grid_time
                    pVnext_given_V_L = zeros([n_grid_voltage n_grid_voltage 2]);
                    pVnext_given_V_L(index_reset,:,2) = 1;
                    % We need to calculte this for each time point since the current changes
                    % each time
                    for i_v = 1:n_grid_voltage
                        v_noise = v_grid-v_grid(i_v)-((E_L-v_grid(i_v))*g + I_stimuli(i_t)*k)*dt;
                        relevant_index = (v_noise > (stoc_mu-3*stoc_sigma)) & (v_noise < (stoc_mu+3*stoc_sigma));
                        pVnext_given_V_L(relevant_index,i_v,1) = normpdf(v_noise(relevant_index),stoc_mu,stoc_sigma)*(v_grid(2)-v_grid(1));
                    end
                    %pVnext_given_VL(v,vp,ip)
                    %pL_given_V(i,v)
                    %pVL_given_I(i_t-1,vp,ip)
                    for i_v = 1:n_grid_voltage
                        temp_II_and_III = pVnext_given_V_L(i_v,:,:).*pVL_given_I(i_t-1,:,:);
                        for i = 1:2
                            pVL_given_I(i_t, i_v, i)=pL_given_V(i,i_v)*sum(sum(temp_II_and_III));
                        end
                    end
                end
                M_intensity{i_trial,i_cell} = sum(pVL_given_I(:,:,2),2);
            end
        end
    end
    fprintf('%d\n', i_cell);
end
