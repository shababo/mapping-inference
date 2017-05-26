%rng(12242,'twister');
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));

i_rep = 100;
factor=20;
num_sample = 50;

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

v_reset_known=-4000;
params_sim.V_th= 15;
params_sim.V_reset = -4000;
params_sim.g =  l23_cells_for_sim(i_template).g;
params_sim.gain =  l23_cells_for_sim(i_template).optical_gain*factor;
funcs.invlink = @invlink_test;%@(resp) log(1 + exp(resp));%@(x) exp(x);
stoc_params.mu=0;
stoc_params.sigma = 1e-8;

%stimuli_seq = num_sample*(rand([num_sample 1])+0.5);
stimuli_seq = 500*rand([num_sample 1])/factor;

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
     [V_vect, spikes]  = lif_glm_sim_v2(stim,params_sim,funcs);
    responses(i_trial,:)=spikes;
    sum(spikes)
end

%
in_params.g =   l23_cells_for_sim(i_template).g;

% LIF-GLM fits
%-------------------------------------%
[stats_conv] = fit_lifglm_v2(responses, stims,in_params,v_reset_known);

% Output:
stats_conv.beta(1)
params_sim.gain
(stats_conv.beta-params_sim.gain)/params_sim.gain
