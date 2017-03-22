%% Fit the LIF-GLM model with soft-assignments
% Updates: allow for more than two spikes, and allow for offsef in the GLM
% fit

%%%%%location dummy variable ind_allLoc

%% Load the soft assignments
n_trials_this_cell = length(trials_precell{i_cell});
spTrain_lifglm = zeros(length(t_vect),n_trials_this_cell);
this_stimuli = zeros(n_trials_this_cell,1);


for i_trial = 1:n_trials_this_cell
    ind_trial =  times_precell{i_cell}(2,:)==trials_precell{i_cell}(i_trial);
    this_assignments = events_precell{i_cell}(2,ind_trial);
    this_time = times_precell{i_cell}(1,ind_trial);
    % Reformate it into discrete vector:
    spTrain_lifglm(round(this_time ,0)+1,i_trial) = this_assignments.*( this_assignments>0.01);
    this_stimuli(i_trial) = stimuli_size_local(trials_precell{i_cell}(i_trial), i_cell);
end

%%
% Do nothing if there are not enough events:
if sum(sum(spTrain_lifglm)) > 10
    
    
    marks = this_stimuli;
    I_eg=repmat(I_stimuli,1,n_trials_this_cell);
    
    g=0.1;
    [cov1,cov2,cov3]=gconv_v3(I_eg,spTrain_lifglm,g,marks); %temporally convolve paramters with g upto spike time
    
    design=[-ones(size(cov1(:))) cov1(:) cov2(:) cov3(:)];
    %design=[ones(size(cov1(:))) cov2(:) cov3(:)];
    Y=spTrain_lifglm(:);
    
    
    num_params=size(design,2)-1;
    init_vals=zeros(1,num_params);
    
    x0 = init_vals; x0(1)=-5;x0(2)=-20;x0(3)=-50;
    
    obj_fun = @(x) logL_LIFglmnn_offset( x, design, Y, link_function);
    options = optimset('Algorithm','interior-point','Display','iter','GradObj','on','Diagnostics','on','UseParallel','always'); %
    
    ub = Inf*ones(size(x0));
    lb = -Inf*ones(size(x0));
    lb(3:num_params) = -1e-6;
    ub(3:num_params) = 1e6;
    
    betahat_conv_allLoc = fmincon(obj_fun,x0,[],[],[],[],lb,ub,[],options); % V_thres, V_resting, V_reset
    
    
    % otherwise do nothing
    % Use the updated parameters to estimate the intensity
    
    
    V_th = betahat_conv_allLoc(1);
    V_reset = betahat_conv_allLoc(3);
    E_L=betahat_conv_allLoc(2); %resting membrane potential [mV]
    
    V_max = V_th+5;
    V_min= V_reset-5;
    v_grid =  (V_max -V_min)*(0: (n_grid_voltage-1) )/(n_grid_voltage-1) +V_min;
    
    [~, index_reset] = min(abs(v_grid-V_reset));
    [~, index_rest] = min(abs(v_grid-E_L));
    
    % Obtain stimuli information on this cell
    if sum(stimuli_size_local(:,i_cell)> k_minimum ) > 5
        
        stim_seq =stimuli_size_local(stimuli_size_local(:,i_cell)> k_minimum,i_cell);
        gap_stimuli=range(stim_seq)/n_stimuli_grid;
        stimuli_bins{i_cell} = min(stim_seq):gap_stimuli:max(stim_seq);
        mid_points = (stimuli_bins{i_cell}(2:end)+stimuli_bins{i_cell}(1:(n_stimuli_grid)))/2;
        
        for i_stimuli = 1:n_stimuli_grid
            k_temp = mid_points(i_stimuli);
            pL_given_V = zeros([2 n_grid_voltage]);
            
            if link_function == 1 % logit link
                pL_given_V(2,:) = min(1,exp(t_factor*(v_grid - V_th))/(1+exp(t_factor*(v_grid - V_th) )));
                
            elseif link_function == 2 % log link
                pL_given_V(2,:) = min(1, exp(t_factor*(v_grid - V_th) ));
                
            end
            
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
        for i_trial = 1:n_trial
            
            % Initialize the storage
            % Simulated traces for marginal intensity
            k = stimuli_size_local(i_trial, i_cell);
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
            k = stimuli_size_local(i_trial, i_cell);
            if k > k_minimum
                pL_given_V = zeros([2 n_grid_voltage]);
                
                
                if link_function == 1 % logit link
                    pL_given_V(2,:) = min(1,exp(t_factor*(v_grid - V_th))/(1+exp(t_factor*(v_grid - V_th) )));
                    
                elseif link_function == 2 % log link
                    pL_given_V(2,:) = min(1, exp(t_factor*(v_grid - V_th) ));
                    
                end
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
    fprintf('%d updated\n', i_cell);
    
else
    
    fprintf('%d not enough events\n', i_cell);
    
    
end