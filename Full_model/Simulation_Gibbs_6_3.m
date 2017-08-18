addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Load the data set 
load('./Environments/6_3_s2c3_mpp_and_stim_data.mat')
%target_locations =target_locs;
load('./Environments/6_3_cell_locs_reduced_and_fixed.mat')
cell_locations = cell_locs;
n_cell=size(cell_locations,1);
target_locations=cell_locations;
%%
% n_cell = 5;
% cell_locations=cell_locations(1:n_cell,:);
% target_locations=cell_locations;
%% Load the current template 
load('./Environments/chrome-template-3ms.mat');
downsamp=1;max_time=300;
power_level = [50 75 100];
num_power_level=length(power_level);
current_template=template(1:downsamp:max_time);
t_vect= 1:1:max_time;
%% Load the shape template 
load('./Environments/l23_template_cell.mat');

load('./Environments/l23_cells_for_sim.mat');
%l23_average_shape
temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;
shape_template=l23_average_shape;
%%
% Calculate the size of stimuli
cell_params.locations =  cell_locations;
cell_params.shape_gain = ones(n_cell,1);
cell_template = struct();
cell_template.shape= shape_template;
% [pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
%     cell_template,target_locations);
[pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
    cell_template,target_locations);
%% Delay parameters
delay_params.type=2; %1: normal; 2: gamma
delay_params.mean=60;
delay_params.std=15;
delay_params.delayed=true;
delay_params.n_grid=200;
%%
locations_trials = repmat(1:n_cell,[1 20])';
powers_trials = 100*ones(length(locations_trials),1);
n_trial = size(locations_trials,1);
stimuli_size=zeros(n_trial,n_cell);
for l = 1:n_trial
%     for m = 1:size(mpp(l).locations,2)
%         if isnan(mpp(l).locations(m))
%         else
            stimuli_size(l,:) = stimuli_size(l,:)+...
                (pi_target(:,locations_trials(l)).*powers_trials(l))';
%         end
%     end
end
%% Randomly assign values and generate data
rng(12242,'twister');
background_rate=1e-4;
v_th_known=15*ones([n_cell,1]);
v_reset_known=-1e4*ones([n_cell,1]);
gain_truth = 0.01+0.01*(rand([n_cell 1]));
g_truth = 0.02*ones([n_cell,1]);
gamma_truth = (0.3+0.7*rand([n_cell 1])).*(rand([n_cell 1])>0.8);
% gamma_truth = (ones([n_cell 1])).*(rand([n_cell 1])>0.8);
% gamma_truth = (ones([n_cell 1]));

funcs.invlink = @invlink_sig;%@(resp) log(1 + exp(resp));%@(x) exp(x);
stim_threshold=10;
time_max =300;
%% Simulate data:
mpp=struct();
mu_bg = 1/background_rate;
for i_trial = 1:n_trial
    mpp(i_trial).times=[];
    mpp(i_trial).locations=locations_trials(i_trial,:);
    mpp(i_trial).power=powers_trials(i_trial,:);
    
    for i_cell = 1:n_cell
        k=stimuli_size(i_trial,i_cell);
        if k > stim_threshold
        params_sim.V_th= v_th_known(i_cell);
        params_sim.V_reset = v_reset_known(i_cell);
        params_sim.g = g_truth(i_cell);
        params_sim.gain =  gain_truth(i_cell);
        
        stim = current_template*k;
        [V_vect, spikes]  =lif_glm_sim_v2(stim,params_sim,funcs);
%         plot(V_vect)
        if sum(spikes)>0
            if delay_params.type==2
                shape=(delay_params.mean^2)/(delay_params.std^2);
                scale = delay_params.mean/shape;
                delay_vec=round(gamrnd(shape,scale,[sum(spikes) 1]));
            else
                delay_vec=zeros([sum(spikes) 1]);
            end
            spikes_delay =find(spikes)+delay_vec';
            for i = 1:sum(spikes)
                if rand(1) < gamma_truth(i_cell)
                % censoring at time max:
                if spikes_delay(i)<time_max
                    mpp(i_trial).times=[mpp(i_trial).times spikes_delay(i)];
                end
                end
                
            end
        end
        end
    end
    % add background event:
    R = exprnd(mu_bg);
    while R < time_max
         mpp(i_trial).times=[mpp(i_trial).times max(1,round(R))];
        
        R = R+exprnd(mu_bg);
    end
     fprintf('Trial %d simulated;\n', i_trial);
end
%%
length([mpp.times])/n_trial

%%
stim_threshold =10;g=0.02;background_rate = 1e-4;v_th_known=15*ones(n_cell,1);
linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};

gain_grid=0.001*[8:30];gain_prior=normpdf(gain_grid,0.015,0.005);
gamma_grid= 0.1*[0:10]; 
gamma_prior=gamma_grid;gamma_prior(1)=0.7;
gamma_prior(2:end)= (1- gamma_prior(1))/(length(gamma_grid)-1);


% gain_initial= 0.015*ones(n_cell,1);
% gamma_initial = 0.5*ones(n_cell,1);
n_gibbs_sample=200;
gain_initial= [];
gamma_initial = [];

% params.invlink = @invlink_sig;
% params.dlink = @derlink_sig;
% params.link = @link_sig;
% params.dinvlink = @derinvlink_sig;
% linkfunc = {params.link, params.dlink, params.invlink,params.dinvlink};
%% run the Gibbs sampler
[gain_samples, gamma_samples] = Gibbs_first_spike(mpp, ...
    target_locations, cell_locations,...
    current_template, shape_template,delay_params,linkfunc,...
    stim_threshold, g, background_rate,v_th_known,...
    gain_grid, gain_prior, gamma_grid,gamma_prior,...
    gamma_initial,gain_initial,n_gibbs_sample);
%%
figure(1)
plot(gamma_truth+normrnd(0,0.01,[n_cell 1]),mean(gamma_samples,2),'.','MarkerSize',20)
xlim([0,1]);
ylim([0,1]);
figure(2)
plot(gain_truth(gamma_truth>0),mean(gain_samples(gamma_truth>0,:),2),'.','MarkerSize',20)
xlim([0.005,0.025]);
ylim([0.005,0.025]);
