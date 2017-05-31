%%
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));

%%
rng(12242,'twister');
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
i_template = 3;

%
params_sim.V_th= 15;
params_sim.V_reset = -4000;
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
current_template=template(1:downsamp:400);
I_e_vect=current_template;
%
responses = zeros(num_sample, length(I_e_vect));
stims = zeros(num_sample, length(I_e_vect));

for i_trial = 1:num_sample
    k=stimuli_seq(i_trial);
    stim = I_e_vect*k;
    stims(i_trial,:) = stim;
     [V_vect, spikes]  = lif_glm_sim_v2(stim,params_sim,funcs);
    responses(i_trial,:)=spikes;
    sum(spikes)
end
%%
avg_spikes=mean(responses,1);
%plot(avg_spikes)
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

n_grid_voltage =1000;
dt=1;
t_grid = dt:dt:length(I_e_vect);
n_grid_time = length(t_grid);
I_stimuli= I_e_vect;
t_factor=1;
k_minimum=0.001;
cell_params=params_sim;
cell_params.gain_sd=0.01;

delay_params.mean=0;
delay_params.std=0;
sd_range=1;
M_intensity_func = Intensity_v5(stimuli_seq, 5,n_grid_voltage,...
        t_grid,t_factor,k_minimum,...
         cell_params, funcs,...
        I_stimuli,sd_range,delay_params);
%%
%
plot(t_grid,M_intensity_func{1,1},'col','r')
hold on;
plot(avg_spikes,'col','b')
hold on;
%plot(t_grid,M_intensity_func{1,1},'col','k')

hold off;
