% Exploit the fact that V_th, V_reset, and g are set to be the same:
function [M_intensity] = Intensity_v4(stimuli_size, n_stimuli_grid,n_grid_voltage,...
        t_grid,t_factor,k_minimum,...
         cell_params, funcs,...
        I_stimuli, stoc_mu, stoc_sigma)

n_trial_temp = size(stimuli_size,1);
stimuli_bins = cell(size(stimuli_size,2),1);
M_grid_intensity = cell(n_stimuli_grid,1);

% marginal expectation
M_intensity=cell(size(stimuli_size));
n_cell_local = size(stimuli_size,2);

n_grid_time = length(t_grid);
        

% Prepare data:
dt = t_grid(2)-t_grid(1);
t_grid_upp = t_grid+dt/2;
t_grid_low = t_grid-dt/2;

V_th = cell_params.V_th(1);
V_reset = cell_params.V_reset(1);
E_L=0;
g=cell_params.g(1);
V_max = V_th+5;
V_min= V_reset-5;
% We have to use a special grid for now since V_reset is -1e3
v_grid =  [linspace(V_min,-30,n_grid_voltage/2) linspace(-25,V_max,n_grid_voltage/2)];
v_gap = v_grid(2:end)-v_grid(1:n_grid_voltage-1);
v_gap = [v_gap v_gap(end)];
[~, index_reset] = min(abs(v_grid-V_reset));
[~, index_rest] = min(abs(v_grid-E_L));

    stim_cell = stimuli_size.* (ones(size(stimuli_size,1),1)*cell_params.gain');
    
    stim_seq = reshape(stim_cell,[size(stimuli_size,1)*size(stimuli_size,2) 1] );
    stim_seq =stim_cell( stim_seq> k_minimum);
    % Obtain stimuli information on this cell
    if sum(stim_seq> k_minimum ) > 5 & length(unique(stim_seq)) >n_stimuli_grid
        gap_stimuli=range(stim_seq)/n_stimuli_grid;
        stimuli_bins= min(stim_seq):gap_stimuli:max(stim_seq);
        mid_points = (stimuli_bins(2:end)+stimuli_bins(1:(n_stimuli_grid)))/2;
    else
        mid_points =unique(stim_seq);
    end
        
        pL_given_V = zeros([2 n_grid_voltage]);
        pVL_given_I = zeros([n_grid_time n_grid_voltage 2]);
        pVnext_given_V_L = zeros([n_grid_voltage n_grid_voltage 2]);
        
        n_stimuli_temp = length(mid_points);
        for i_stimuli = 1:n_stimuli_temp
            k_temp = mid_points(i_stimuli);
            pL_given_V(:,:)=0;
            pVL_given_I(:,:,:)=0;
            pL_given_V(2,:) =  min(1,t_factor*funcs.invlink(v_grid-V_th));
            
            pL_given_V(1,:) = 1-pL_given_V(2,:);
            pVL_given_I(1,index_rest,1)=1;
            
            normconstant = 1/(sqrt(2*pi)*stoc_sigma);
            for i_t = 2:n_grid_time
                pVnext_given_V_L(:,:,:)=0;
                pVnext_given_V_L(index_reset,:,2) = 1;
                % We need to calculte this for each time point since the current changes
                % each time
                for i_v = 1:n_grid_voltage
                    v_noise = v_grid-v_grid(i_v)-((E_L-v_grid(i_v))*g +...
                        I_stimuli(i_t)*k_temp)*dt;
                    relevant_index = (v_noise > (stoc_mu-3*stoc_sigma)) & (v_noise < (stoc_mu+3*stoc_sigma));
                    %relevant_index = 1:length(v_noise);
                    % faster than the normpdf()..
                    pVnext_given_V_L(relevant_index ,i_v,1) = ...
                      exp(-(v_noise(relevant_index)-stoc_mu).^2/(2*stoc_sigma^2)).*v_gap(relevant_index)*normconstant;
                end
                for i_v = 1:n_grid_voltage
                    temp_II_and_III = pVnext_given_V_L(i_v,:,1)*pVL_given_I(i_t-1,:,1)';
                    pVL_given_I(i_t, i_v, :)=pL_given_V(:,i_v)*temp_II_and_III;
                end
                % Fixing the reset probability:
                pVL_given_I(i_t, index_reset, 1)=pVL_given_I(i_t, index_reset, 1)+...
                        pL_given_V(1,index_reset)*(pVnext_given_V_L(index_reset,:,2)*pVL_given_I(i_t-1,:,2)');
            end
            M_grid_intensity{i_stimuli} = sum(pVL_given_I(:,:,2),2);
            
        fprintf('%d grid \n',i_stimuli);
        end
        
        for i_cell = 1:n_cell_local
            
            
            for i_trial = 1:n_trial_temp
                % Initialize the storage
                % Simulated traces for marginal intensity
                k = stimuli_size(i_trial, i_cell)*cell_params.gain(i_cell);
                M_intensity{i_trial,i_cell} = zeros([length(t_grid) 1]);
                if k > k_minimum
                    [~, ind_seq] = min(abs(k-mid_points));
                    M_intensity{i_trial,i_cell} = M_grid_intensity{ind_seq};
                end
            end
            fprintf('%d\n', i_cell);
        end
