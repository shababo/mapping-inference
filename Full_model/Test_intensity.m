%% Testing the intensity estimation 
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%%
rng(12242,'twister');
num_sample = 5000;
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
time_max= 300;
current_template=template(1:downsamp:time_max);
I_e_vect=current_template;
%

% delay_params.mean=30;
% delay_params.std=5;

%delay = normrnd(delay_params.mean,delay_params.std);
%delay = 40;


n_grid_voltage =1000;
stimulus_threshold=-50;
dt=1;
n_grid_time = length(I_e_vect);
I_stimuli= I_e_vect;
t_factor=1;
k_minimum=0.001;
cell_params=params_sim;
cell_params.gain_sd=0.01;

sd_range=1;

%delay_params.mean=round(delay);
%delay_params.std=0.01;
n_delay_grid = 200;

stimulus_threshold=0.1;
n_stimuli_grid=0;
gap_stimuli=5;

stimuli_size=stimuli_seq;

outputM=true;

first_only=true;
V_threshold = -50;
%%
for i_outer = 1:3
    responses = zeros(num_sample, length(I_e_vect));
    stims = zeros(num_sample, length(I_e_vect));
    
    if i_outer==1
        delay_params.type=1;
        delay_params.mean=0;
        delay_params.std=0.1;
    elseif i_outer==2
        delay_params.type=0;
        delay_params.shape=0.8*30*30/25;
        delay_params.scale=25/30;
    else
        delay_params.type=1;
        delay_params.mean=30;
        delay_params.std=5;
    end
    for i_trial = 1:num_sample
        k=stimuli_seq(i_trial);
        stim = I_e_vect*k;
        stims(i_trial,:) = stim;
        [V_vect, spikes]  = lif_glm_sim_v2(stim,params_sim,funcs);
        
        % shift the spikes give the delay
        
        if i_outer==1
            delay=0;
        elseif i_outer==2
            delay = round(gamrnd(delay_params.shape,delay_params.scale,[sum(spikes) 1]));
        elseif i_outer==3
            
            delay = round(normrnd(delay_params.mean,delay_params.std,[sum(spikes) 1]));
        end
        
        
        spikes_delayed=find(spikes)+round(delay)';
          spikes_delayed=spikes_delayed( spikes_delayed<time_max);
        responses(i_trial,spikes_delayed)=1;
        %sum(responses(i_trial,:))
    end
    
    mpp=struct();
    for i=1:num_sample
        mpp(i).times=find(responses(i,:));
    end
    avg_spikes{i_outer}=mean(responses,1);
    
    
    
    [A,B,M_intensity_func]=Intensity_v7(outputM,... % whether to output M intensity
        stimuli_size, mpp,I_stimuli,... % data from exp
        cell_params, delay_params,... % estimated parameters
        funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
        n_stimuli_grid,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
        V_threshold,stimulus_threshold,first_only);
    estimated_intensities{i_outer}=M_intensity_func{1,1};
end
%%
%
figure(1)
colors=['r' 'g' 'b'];
t_grid = 1:length(I_e_vect);
for i_outer = 1:3
    plot(t_grid,estimated_intensities{i_outer},'col',colors(i_outer),'LineStyle','--','LineWidth',1)
    hold on;
    plot(t_grid,avg_spikes{i_outer},'col',colors(i_outer),'LineWidth',1)
    hold on;
end
xlabel('Time (1/20 ms)')
ylabel('Intensities')
hold off;
% saveas(1,'./Intensity_check.png');