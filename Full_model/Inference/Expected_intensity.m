%% Simulate multiple presynaptic processes
% And calculate the expected intensity
%% Parameters:
%
if exact_crossing == 1
    simulated_events = cell(n_trial, n_cell);
    E_intensity=cell(n_trial, n_cell);
    
    
    n_sim = 50;
    simulated_events_index = zeros( [n_sim length(t_vect)]);
    
    for i_cell = 1:n_cell
        
        V_th = all_V_th(i_cell);
        V_reset = all_V_reset(i_cell);
        E_L=all_E_L(i_cell); %resting membrane potential [mV]
        num_I_Stim=1; % Simulation scheme
        
        
        % Simulated traces for marginal intensity
        %
        for i_trial = 1:n_trial
            simulated_events_index = zeros( [n_sim length(t_vect)]);
            simulated_events{i_trial, i_cell} = [];
            t_vect=0:data_params.dt:data_params.T;
            k = stimuli_size(i_trial, i_cell);
            if k > 0.00001
                for i_sim = 1:n_sim
                    
                    V_vect=zeros(1,length(t_vect));
                    
                    %%%DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS
                    
                    %spTrain=zeros(t_end,length(I_Stim_vect));% The output spike train
                    i=1; %index denoting which element of V is being assigned
                    V_vect(i)=E_L; %first element of V, i.e. value of V at t=0
                    %%%%chi-sq shape current
                    I_e_vect=[0;I_e(:,num_I_Stim)];
                    for ti=data_params.dt:data_params.dt:data_params.T %loop through values of t in steps of df ms
                        V_vect(i+1) = V_vect(i) + ((E_L-V_vect(i))*g + I_e_vect(i)*k)*data_params.dt + ...
                            sqrt(data_params.dt)*normrnd(stoc_mu,stoc_sigma);
                        %if statement below says what to do if voltage crosses threshold
                        if (V_vect(i+1)>V_th) %cell spiked
                            V_vect(i+1)=V_reset; %set voltage back to V_reset
                            %simulated_events{i_trial, i_cell} = [simulated_events{i_trial, i_cell} i];
                            simulated_events_index(i_sim, i)=1;
                            %t
                        end
                        i=i+1; %add 1 to index,corresponding to moving forward 1 time step
                    end
                    
                end
                
            end
            % Empirical intensity
            E_intensity{i_trial, i_cell} = zeros(length(t_vect),1);
            E_intensity{i_trial, i_cell}=sum(simulated_events_index,1)/n_sim;
        end
        fprintf('%d\n', i_cell);
    end
else
    % marginal expectation
    M_intensity=cell(n_trial, n_cell);
    for i_cell = 1:n_cell
        % Initialize the storage
        n_grid_voltage = 200;
        n_grid_time = length(t_vect);
        t_grid=t_vect;
        t_factor = 1; % dt = 1/20 ms.
        V_max = 0;
        v_grid =  (V_max - all_V_reset(i_cell))*(0: (n_grid_voltage-1) )/(n_grid_voltage-1) +all_V_reset(i_cell);
        [~, index_reset] = min(abs(v_grid-all_V_reset(i_cell)));
        [~, index_rest] = min(abs(v_grid-all_E_L(i_cell)));
        
        V_th = all_V_th(i_cell);
        V_reset = all_V_reset(i_cell);
        E_L=all_E_L(i_cell); %resting membrane potential [mV]
        num_I_Stim=1; % Simulation scheme
        
        
        % Simulated traces for marginal intensity
        %
        for i_trial = 1:n_trial
            % Prepare data:
            t_grid_upp = t_grid+data_params.dt/2;
            t_grid_low = t_grid-data_params.dt/2;
            pL_given_V = zeros([2 n_grid_voltage]);
            pL_given_V(2,:) = min(1,t_factor*max(v_grid - all_V_th(i_cell),0));
            pL_given_V(1,:) = 1-pL_given_V(2,:);
            
            pVL_given_I = zeros([n_grid_time n_grid_voltage 2]);
            pVL_given_I(1,index_rest,1)=1;
            
            k = stimuli_size(i_trial, i_cell);
            if k > 0.00001
                for i_t = 2:n_grid_time
                    pVnext_given_V_L = zeros([n_grid_voltage n_grid_voltage 2]);
                    pVnext_given_V_L(index_reset,:,2) = 1;
                    % We need to calculte this for each time point since the current changes
                    % each time
                    for i_v = 1:n_grid_voltage
                        v_noise = v_grid-v_grid(i_v)-((E_L-v_grid(i_v))*g + I_e_vect(i_t)*k)*data_params.dt;
                        pVnext_given_V_L(:,i_v,1) = normpdf(v_noise,stoc_mu,stoc_sigma)*(v_grid(2)-v_grid(1));
                        % Might be possible to reduce the comp cost by considering only a few
                        % neighbours
                        %sum(pVnext_given_V_L(:,i_v,1))
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
                
                
            end
            M_intensity{i_trial, i_cell} = sum(pVL_given_I(:,:,2),2);
        end
        fprintf('%d\n', i_cell);
    end
end

