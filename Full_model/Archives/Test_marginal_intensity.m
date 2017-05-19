rng(12242,'twister');
i_cell = 140;
i_cell_local = 29;

stimuli_size_copy = stimuli_size_local(:,i_cell_local);

%max(max( abs(stimuli_size(:,local_index) - stimuli_size_local)))
%k_basic*unifrnd(0,1,n_trial,1);

presynaptic_events = cell(n_trial, 1);
presynaptic_amplitudes = cell(n_trial, 1);
voltage_traces = cell(n_trial, 1);
% One cell:
i_cell = 140;
V_th = all_V_th(i_cell);
V_reset = all_V_reset(i_cell);
E_L=all_E_L(i_cell); %resting membrane potential [mV]
num_I_Stim=1;
for i_trial = 1:n_trial
    presynaptic_events{i_trial, 1} = [];
    presynaptic_amplitudes{i_trial,1} = [];
    k = stimuli_size(i_trial,i_cell);
    if k > k_minimum
        %%%DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS
        t_vect=0:data_params.dt:data_params.T;
        V_vect=zeros(1,length(t_vect));
        
        %spTrain=zeros(t_end,length(I_Stim_vect));% The output spike train
        i=1; %index denoting which element of V is being assigned
        V_vect(i)=E_L; %first element of V, i.e. value of V at t=0
        %%%%chi-sq shape current
        I_e_vect=[0;I_e(:,num_I_Stim)];
        for ti=data_params.dt:data_params.dt:data_params.T %loop through values of t in steps of df ms
            V_vect(i+1) = V_vect(i) + ((E_L-V_vect(i))*g + I_e_vect(i)*k)*data_params.dt + ...
                sqrt(data_params.dt)*normrnd(stoc_mu,stoc_sigma);
            
                prob_spike = min(1,max(0,V_vect(i+1) - V_th)*data_params.dt);
                spike_indicator = binornd(1,prob_spike);
                if spike_indicator %cell spiked
                    V_vect(i+1)=V_reset; %set voltage back to V_reset
                    presynaptic_events{i_trial, 1} = [presynaptic_events{i_trial, 1} t_vect(i)];
                    presynaptic_amplitudes{i_trial, 1} = [presynaptic_amplitudes{i_trial, 1} normrnd(all_amplitudes(i_cell),evoked_params.sigma_a)];
                end
            i=i+1; %add 1 to index,corresponding to moving forward 1 time step
        end
        voltage_traces{i_trial,1} = V_vect;
    end 
end
fprintf('%d\n', i_cell);
%%
all_events = [presynaptic_events{:,1}];
size(all_events)
%%
load('./Environments/current_template.mat'); %Contains the vector norm_average_current
I_eg_100mV=abs(norm_average_current(1:20:20*75).*100);
I_eg_50mV=abs(norm_average_current(1:20:20*75).*50);
I_eg_25mV=abs(norm_average_current(1:20:20*75).*25);
I_e=[repmat(I_eg_100mV',1,5) repmat(I_eg_50mV',1,5) repmat(I_eg_25mV',1,5)];
num_I_Stim=1;
I_e_vect=[0;I_e(:,num_I_Stim)];
I_stimuli = I_e_vect;
% Stochastic components of voltages 
stoc_mu=0;stoc_sigma=0.5;
g=0.1; %membrane time constant [ms]

T=75;
dt=1;
t_vect=0:dt:T;

V_thresholds = local_V_th;
V_resets = local_V_reset;
E_Ls = local_E_L;

n_stimuli_grid=40;
k_basic = 0.04;
exact_crossing = 0;
%% Estimate the marginal intensity:
% Calculate ten possible trajectories for each cell
% at 1, 0.9, 0.8, ..., 0.1 amount of stimulation
% We can even do this using the actual stimulation each cell received
stimuli_bins = cell(1,1);

M_grid_intensity = cell(n_cell_local,n_stimuli_grid);
% marginal expectation
M_intensity=cell(n_trial,1);
n_grid_voltage = 400;
n_grid_time = length(t_vect);
t_grid=t_vect;

stoc_sigma= 0.5;
t_factor = 1;

% Prepare data:
t_grid_upp = t_grid+dt/2;
t_grid_low = t_grid-dt/2;

V_th = V_thresholds(i_cell);
V_reset = V_resets(i_cell);
E_L=E_Ls(i_cell); %resting membrane potential [mV]

%V_max = V_th+max(1/t_factor,3*stoc_sigma);
%V_max= 0;
V_max = V_th+5;
V_min= V_reset-5;
v_grid =  (V_max -V_min)*(0: (n_grid_voltage-1) )/(n_grid_voltage-1) +V_min;

[~, index_reset] = min(abs(v_grid-V_reset));
[~, index_rest] = min(abs(v_grid-E_L));

% Obtain stimuli information on this cell
%sum(stimuli_size(:)> k_minimum)


stim_seq =stimuli_size_copy(stimuli_size_copy> k_minimum);
gap_stimuli=range(stim_seq)/n_stimuli_grid;
stimuli_bins{i_cell} = min(stim_seq):gap_stimuli:max(stim_seq);
mid_points = (stimuli_bins{i_cell}(2:end)+stimuli_bins{i_cell}(1:(n_stimuli_grid)))/2;

%%

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
    
    M_grid_intensity{1, i_stimuli} = sum(pVL_given_I(:,:,2),2);
end

%% Using simulation
n_sim= 200;
M_grid_intensity = cell(1,n_stimuli_grid);
% marginal expectation
M_intensity=cell(n_trial,1);

for i_stimuli = 1:n_stimuli_grid
    k_temp = mid_points(i_stimuli);
    t_vect=0:data_params.dt:data_params.T;
    simulated_events_index= zeros( [n_sim length(t_vect)]);
    for i_sim = 1:n_sim
        V_vect=zeros(1,length(t_vect));
        %spTrain=zeros(t_end,length(I_Stim_vect));% The output spike train
        i=1; %index denoting which element of V is being assigned
        V_vect(i)=E_L; %first element of V, i.e. value of V at t=0
        for ti=data_params.dt:data_params.dt:data_params.T %loop through values of t in steps of df ms
            V_vect(i+1) = V_vect(i) + ((E_L-V_vect(i))*g + I_e_vect(i)*k_temp)*data_params.dt + ...
                sqrt(data_params.dt)*normrnd(stoc_mu,stoc_sigma);
            prob_spike = min(1,max(0,V_vect(i+1) - V_th)*data_params.dt);
            spike_indicator = binornd(1,prob_spike);
            if spike_indicator %cell spiked
                V_vect(i+1)=V_reset; %set voltage back to V_reset
                simulated_events_index(i_sim, i) = 1;
                %t
            end
            i=i+1; %add 1 to index,corresponding to moving forward 1 time step
        end
    end
    M_grid_intensity{1, i_stimuli} =  mean(simulated_events_index,1);
end

%%

for i_trial = 1:n_trial
    
    % Initialize the storage
    % Simulated traces for marginal intensity
    k = stimuli_size_copy(i_trial);
    if k > k_minimum
        [~, ind_seq] = min(abs(k-mid_points));
        M_intensity{i_trial,1} = M_grid_intensity{1, ind_seq};
    else
        M_intensity{i_trial,1} = zeros([length(t_vect) 1]);
    end
end
%%
estimated_intensity = M_intensity;
expected_all = zeros(i_trial,1);
for i_trial = 1:n_trial
    expected_all(i_trial,1)=sum(estimated_intensity{i_trial, 1} );
end

%%
size([presynaptic_events{:,1}])/sum(expected_all)

%% 
expected_sim = expected_all

%%
sum(expected_all)-sum(expected_sim)