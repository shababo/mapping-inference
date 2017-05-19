%%
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));

%rng(12242,'twister');
factor=20;
num_sample = 500;

load('./Environments/l23_cells_for_sim.mat');

num_types_cell = length(l23_cells_for_sim);
% normalized the cell shapes
for i = 1:num_types_cell
    temp=l23_cells_for_sim(i).shape;
    temp_max = max(max(max(temp)));
    l23_cells_for_sim(i).shape = temp/temp_max;
end
% 
i_template = 1;

%
params_sim.V_th= 15;
params_sim.V_reset = -2000;
params_sim.g =  l23_cells_for_sim(i_template).g;
params_sim.gain =  l23_cells_for_sim(i_template).optical_gain;
funcs.invlink = @invlink_test;%@(resp) log(1 + exp(resp));%@(x) exp(x);
stoc_params.mu=0;
stoc_params.sigma = 1e-8;

%stimuli_seq = num_sample*(rand([num_sample 1])+0.5);
stimuli_seq = 100*ones(num_sample,1);

%
load('./Environments/chrome-template-3ms.mat');
downsamp=1;
current_template=template(1:downsamp:200);
I_e_vect=current_template;
%
responses = zeros(num_sample, length(I_e_vect));
stims = zeros(num_sample, length(I_e_vect));

for i_trial = 1:num_sample
    k=stimuli_seq(i_trial);
    stim = I_e_vect*k;
    stims(i_trial,:) = stim;
     [V_vect, spikes]  = lif_glm_sim_v2(stim,params_sim,funcs,stoc_params);
    responses(i_trial,:)=spikes;
    sum(spikes)
end
%%
avg_spikes=mean(responses,1);
plot(avg_spikes)
% %
% in_params.g =   l23_cells_for_sim(i_template).g;
% 
% % LIF-GLM fits
% %-------------------------------------%
% [stats_conv] = fit_lifglm_v2(responses, stims,in_params);
% 
% % Output:
% stats_conv.beta(1)
% stats_conv.beta(2)
% params_sim.gain
% (stats_conv.beta(2)-params_sim.gain)/params_sim.gain

%%
n_grid_voltage =2000;
dt=1;
t_grid = dt:dt:length(I_e_vect);
n_grid_time = length(t_grid);
I_stimuli= I_e_vect;

% Prepare data:
dt = t_grid(2)-t_grid(1);
t_grid_upp = t_grid+dt/2;
t_grid_low = t_grid-dt/2;


V_th =params_sim.V_th;
V_reset = params_sim.V_reset;
E_L=0;
g=params_sim.g;
gain = params_sim.gain;
V_max = V_th+10;
V_min= V_reset-10;



% We have to use a special grid for now since V_reset is -1e3
v_grid =  [linspace(V_min,-30,n_grid_voltage/2) linspace(-25,V_max,n_grid_voltage/2)];
[~, index_reset] = min(abs(v_grid-V_reset));
[~, index_rest] = min(abs(v_grid-E_L));





pL_given_V = zeros([2 n_grid_voltage]);
pVL_given_I = zeros([n_grid_time n_grid_voltage 2]);
pVnext_given_V_L = zeros([n_grid_voltage n_grid_voltage 2]);
        
t_factor=1;


k_temp = stimuli_seq(1)*gain;

pL_given_V(:,:)=0;
pVL_given_I(:,:,:)=0;
pL_given_V(2,:) =  min(1,t_factor*funcs.invlink(v_grid-V_th));

pL_given_V(1,:) = 1-pL_given_V(2,:);
pVL_given_I(1,index_rest,1)=1;
            
for i_t = 2:n_grid_time
    pVnext_given_V_L(:,:,:)=0;
    pVnext_given_V_L(index_reset,:,2) = 1;
    % We need to calculte this for each time point since the current changes
    % each time
    for i_v = 1:n_grid_voltage
        v_noise = v_grid-v_grid(i_v)-((E_L-v_grid(i_v))*g +...
            I_stimuli(i_t)*k_temp)*dt;
        
        [~, relevant_index] = min(abs(v_noise)); 
        %relevant_index = 1:length(v_noise);
        % faster than the normpdf()..
        pVnext_given_V_L(relevant_index ,i_v,1) = 1;
    end
    for i_v = 1:n_grid_voltage
        temp_II_and_III = pVnext_given_V_L(i_v,:,1)*pVL_given_I(i_t-1,:,1)';
        pVL_given_I(i_t, i_v, :)=pL_given_V(:,i_v)*temp_II_and_III;
    end
    % Fixing the reset probability:
    pVL_given_I(i_t, index_reset, 1)=pVL_given_I(i_t, index_reset, 1)+...
        pL_given_V(1,index_reset)*(pVnext_given_V_L(index_reset,:,2)*pVL_given_I(i_t-1,:,2)');
end
M_intensity = sum(pVL_given_I(:,:,2),2);


%
plot(t_grid,M_intensity,'col','r')
hold on;
plot(avg_spikes,'col','b')
hold off;