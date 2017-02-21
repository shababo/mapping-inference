% Gibbs sampler on the full model w/ synpatic failure
% Note:
%   - Sample forward with backward filter algorithm
%   - Assignments

gamma_current = 0.8;
i_trial=1;
% Initialize the storage
n_grid_voltage = 200;
n_grid_time = length(t_vect);
t_grid=t_vect;
t_factor = 1/20; % dt = 1/20 ms.
V_max = 0;
v_grid =  (V_max - all_V_reset(i_cell))*(0: (n_grid_voltage-1) )/(n_grid_voltage-1) +all_V_reset(i_cell);
[~, index_reset] = min(abs(v_grid-all_V_reset(i_cell)));
[~, index_rest] = min(abs(v_grid-all_E_L(i_cell)));


% Prepare data:
t_grid_upp = t_grid+data_params.dt/2;
t_grid_low = t_grid-data_params.dt/2;
d_grid = t_grid;
for i_d = 1:length(t_grid)
    d_grid(i_d) =  sum(assigned_events < t_grid_upp(i_d) & assigned_events > t_grid_low(i_d))>0;
end



%pD_given_V_L(k+2)*pD_given_L(k+1)*pVnext_given_V_L*pL_given_V

pD_given_L = zeros([2 2]);
% First column means L = 0, second column means L=1;
pD_given_L(1,1) = 1;
pD_given_L(2,2) = gamma_current;
pD_given_L(1,2) = 1-gamma_current;
% Second column means L =1

pL_given_V = zeros([2 n_grid_voltage]);
pL_given_V(2,:) = min(1,t_factor*max(v_grid - all_V_th(i_cell),0));
pL_given_V(1,:) = 1-pL_given_V(2,:);

% The recursive algorithm
% The backward probability
pD_given_V_L = zeros([n_grid_time n_grid_voltage 2]);
% Initialization

pD_given_V_L(n_grid_time,:,:)=1/(n_grid_voltage*2);
for i_t_counter = 1:(n_grid_time-1)
    i_t = n_grid_time- i_t_counter;
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
    
    % pD_given_V_L(i_t+1,vp,ip) % I in the integral
    % pL_given_V(ip,vp) % VI in the integral
    
    % pD_given_L(d_grid(i_t)+1, i); %select the row using D(i_t), II
    % pVnext_given_V_L(vp,v,i); % III
    
    
    % This two terms have nothing to do with the value of the current V and
    % L
    temp_I_and_IV = reshape(pD_given_V_L(i_t+1,:,:), [n_grid_voltage 2]).*transpose(pL_given_V(:,:));
    
    temp_I_and_IV_overL = sum(temp_I_and_IV,2);
    
    temp_II_and_III = zeros([n_grid_voltage n_grid_voltage 2]);
    temp_II_and_III(:,:,1) = pVnext_given_V_L(:,:,1)*pD_given_L(d_grid(i_t)+1, 1);
    temp_II_and_III(:,:,2) = pVnext_given_V_L(:,:,2)*pD_given_L(d_grid(i_t)+1, 2);
    
    for i_v = 1:n_grid_voltage
        for i = 1:2
            pD_given_V_L(i_t,i_v,i)  = sum(temp_I_and_IV_overL.*temp_II_and_III(:,i_v,i));
        end
    end
    
    %fprintf('Max prob no fire: %d \n ', max(pD_given_V_L(i_t,:,1)));
    %fprintf('Max prob fire: %d \n', max(pD_given_V_L(i_t,:,2)));
    
end


n_sim_fb = 10;
L_seq = zeros([n_grid_time n_sim_fb]);
V_index = zeros([n_grid_time n_sim_fb]);
V_seq = zeros([n_grid_time n_sim_fb]);

% Draw the first one using prior knowledge:

for i_sim = 1:n_sim_fb
    V_index(1,i_sim) = index_rest; % should be a distribution
    L_seq(1,i_sim)=0;
    V_seq(1,i_sim) = v_grid( index_rest);
    % Draw the rest
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
        end
        
        %pL_given_V(ip,vp)
        %pVnext_given_V_L(vp,v,i)
        %pD_given_V_L(i_t, vp, ip) % the next one
        
        % Calculate the probability
        pVLnext_given_V_L_D = zeros([n_grid_voltage 2]);
        pVLnext_given_V_L_D = transpose(pL_given_V).*reshape(pD_given_V_L(i_t,:,:), [n_grid_voltage 2]).*...
            [pVnext_given_V_L(:,V_index(i_t-1,i_sim),L_seq(i_t-1,i_sim)+1) pVnext_given_V_L(:,V_index(i_t-1,i_sim),L_seq(i_t-1,i_sim)+1)];
        
        % Sample from the derived probability
        
        pL_mar = sum(pVLnext_given_V_L_D)/ sum(sum(pVLnext_given_V_L_D));
        
        L_seq(i_t,i_sim)=binornd(1,pL_mar(2));
        pV_cond = pVLnext_given_V_L_D(:,L_seq(i_t,i_sim)+1);
        pV_cond = pV_cond/sum(pV_cond); % normalize the prob
        V_index(i_t,i_sim)=find(mnrnd(1,pV_cond)); % report the index
        V_seq(i_t,i_sim) = v_grid(V_index(i_t,i_sim));
    end
end





