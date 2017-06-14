% Exploit the fact that V_th, V_reset, and g are set to be the same:
% allow for synaptic delay (distribution or single value)
% Record the number of expected spikes, and the firing rate at the event
% time

% Outputting two kinds of intensities:
%  - Intensities over time 
%  - Intensities convoluted with delay only at the event times 
function [Stimuli_grid, Intensity_grid] = Intensity_v8(... % whether to output M intensity
    stimuli_size, mpp,I_stimuli,... % data from exp
    cell_params,... % estimated parameters
    funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
    n_stimuli_grid,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
    V_threshold,stimulus_threshold,first_only) % thresholds

n_trial = size(stimuli_size,1);
n_cell = size(stimuli_size,2);

M_grid_intensity = cell(n_stimuli_grid,1);

n_grid_time = length(I_stimuli);
dt = 1;
t_grid = 1:dt:n_grid_time;
t_grid_upp = t_grid+dt/2;
t_grid_low = t_grid-dt/2;

V_th = cell_params.V_th(1); % the same threshold
V_reset = cell_params.V_reset(1); % not fitting the reset
E_L=0;

g=cell_params.g(1);
V_max = V_th+5;
V_min= V_reset-5;


stim_cell = stimuli_size.*(ones(size(stimuli_size,1),1)*cell_params.gain');

stim_seq = reshape(stim_cell,[size(stimuli_size,1)*size(stimuli_size,2) 1] );
stim_seq =stim_cell( stim_seq > stimulus_threshold);
% Obtain stimuli information on this cell

if n_stimuli_grid == 0
    n_stimuli_grid =  max(1,round(range(stim_seq)/gap_stimuli));
    fprintf('Stimuli grid change to %d\n', n_stimuli_grid);
else
    gap_stimuli= range(stim_seq)/n_stimuli_grid; 
end

if length(stim_seq) > 5 & length(unique(stim_seq)) >n_stimuli_grid
    gap_stimuli= range(stim_seq)/n_stimuli_grid;
    stimuli_bins= min(stim_seq):gap_stimuli:max(stim_seq);
    mid_points = (stimuli_bins(2:end)+stimuli_bins(1:(n_stimuli_grid)))/2;
else
    mid_points =unique(stim_seq);
end


% Let the gap be the minimum increase in the stimuli input:
if first_only
    v_gap_low = (V_threshold - V_min)/(n_grid_voltage/2);
%    v_gap_upp = (V_max-V_threshold)/(n_grid_voltage/2);

	    v_gap = min(mid_points)*I_stimuli(end);
    v_gap_upp = min(v_gap, (V_max-V_threshold)/(n_grid_voltage/2));

else
    v_gap = min(mid_points)*I_stimuli(end);

    v_gap_low = min(v_gap, (V_threshold - V_min)/(n_grid_voltage/2));
    v_gap_upp = min(v_gap, (V_max-V_threshold)/(n_grid_voltage/2));

end

v_grid =  [V_threshold:v_gap_upp:V_max];
v_grid = [V_min:v_gap_low:V_threshold  v_grid(2:end)];

[~, index_reset] = min(abs(v_grid-V_reset));
[~, index_rest] = min(abs(v_grid-E_L));
if length(v_grid)>n_grid_voltage
    n_grid_voltage =length(v_grid);
     fprintf('Voltage grid change to %d\n', n_grid_voltage);
end
 
v_gap = v_grid(2:end)-v_grid(1:n_grid_voltage-1);
v_gap = [v_gap v_gap(end)];

pL_given_V = zeros([2 n_grid_voltage]);



pL_given_V(2,:) =  min(1,funcs.invlink(v_grid-V_th));
pL_given_V(1,:) = 1-pL_given_V(2,:);
    
parfor i_stimuli = 1:length(mid_points)
    k_temp = mid_points(i_stimuli);
    %reset the values 
    pVL_given_I = zeros([n_grid_time n_grid_voltage 2]);
    pVnext_given_V_L = zeros([n_grid_voltage n_grid_voltage 2]);
    pVL_given_I(:,:,:)=0;
    pVL_given_I(1,index_rest,1)=1;
    for i_t = 2:n_grid_time
        pVnext_given_V_L(:,:,:)=0;
        pVnext_given_V_L(index_reset,:,2) = 1;
        % the transition matrix 
        for i_v = 1:n_grid_voltage
            v_noise = v_grid-...
                (v_grid(i_v)+((E_L-v_grid(i_v))*g +I_stimuli(i_t-1)*k_temp)*dt);
            [~, relevant_index] = min(abs(v_noise));
            pVnext_given_V_L(relevant_index ,i_v,1) = 1;
            if first_only & i_v-1 <index_reset
                 pVnext_given_V_L(relevant_index,i_v,1)=0;
            end
        end
        for i_v = 1:n_grid_voltage
            temp_II_and_III = pVnext_given_V_L(i_v,:,1)*pVL_given_I(i_t-1,:,1)';
            pVL_given_I(i_t, i_v, :)=pL_given_V(:,i_v)*temp_II_and_III;
        end
        % Fixing the reset probability:
        pVL_given_I(i_t, index_reset, 1)=pVL_given_I(i_t, index_reset, 1)+...
            pL_given_V(1,index_reset)*(pVnext_given_V_L(index_reset,:,2)*pVL_given_I(i_t-1,:,2)');
    end
    intensity_temp = sum(pVL_given_I(:,:,2),2);
    M_grid_intensity{i_stimuli} = intensity_temp;
   if mod(i_stimuli,5)==0 
    fprintf('%d grid \n',i_stimuli);
   end
end
Stimuli_grid =mid_points;
Intensity_grid=M_grid_intensity;
