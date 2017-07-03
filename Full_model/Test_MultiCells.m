%% Testing the effect of delay on lif-glm & firing rate
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Read in the cell locations
load('./Environments/6_3_cell_locs_reduced_and_fixed.mat')
local_locations = cell_locs;
%% Draw the cell map
this_cell=36;
figure(1)
temp2 = scatter(local_locations(this_cell,2)+151,-local_locations(this_cell,1)-151,...
    100,'MarkerEdgeColor','g','MarkerFaceColor','b',...
    'MarkerFaceAlpha',0.5);
hold on;
temp1 = scatter(local_locations(:,2)+151,-local_locations(:,1)-151,...
    2);
set(temp1,'MarkerEdgeColor','g','MarkerFaceColor','g');
alpha(temp1,1);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5])
xlim([-20,313]);
ylim([-313,20]);
%% Use cell 133 and its neighbour
dist= zeros(size(local_locations,1),1);
for i= 1:size(local_locations,1)
    dist(i) = sqrt(sum( (local_locations(this_cell,:) -local_locations(i,:)).^2));
end
selected_cells = find(dist<20);
selected_cells
%% A small two cell system:
cell_locations=local_locations(selected_cells,:);
%% Simulation parameters
n_selected_cells = length(selected_cells);
v_th_known=15*ones(length(selected_cells),1);
v_reset_known=-4000*ones(length(selected_cells),1);
gain_initial=0.02*ones(n_selected_cells,1);
gamma=zeros(n_selected_cells,1);
gamma(1)=0.6;
cell_templates = 1:n_selected_cells;

delay_params.type=1;
delay_params.mean=35;
delay_params.std=10;

% background firing rate
background_rate =4e-4;
% length of experiment time (of interest)
time_max= 300;
% number of trials per power levels at each cell
num_per_power=25;
%% Load the cell shape templates for simulation
rng(12242,'twister');
load('./Environments/l23_cells_for_sim.mat');
num_types_cell = length(l23_cells_for_sim);
% normalized the cell shapes
for i = 1:num_types_cell
    temp=l23_cells_for_sim(i).shape;
    temp_max = max(max(max(temp)));
    l23_cells_for_sim(i).shape = temp/temp_max;
end
funcs.invlink = @invlink_sig;%@(resp) log(1 + exp(resp));%@(x) exp(x);
load('./Environments/chrome-template-3ms.mat');
downsamp=1;
current_template=template(1:downsamp:time_max);
%% Other constants:
n_grid_voltage =1000;
voltage_threshold=-50;
dt=1;
n_grid_time = time_max;
t_factor=1;
sd_range=0.1;
n_delay_grid = 200;
stimulus_threshold=1e-4;
n_stimuli_grid=10;
gap_stimuli=0.05;
V_threshold=-50;
first_only=true;
%% Pick the stimulation locations & power
trials_powers = repmat(...
    [50*ones(num_per_power,1); 75*ones(num_per_power,1); 100*ones(num_per_power,1)],...
    [n_selected_cells 1]);
trials_locations = cell_locations(repmat(...
    1:size(cell_locations,1),...
    [3*num_per_power 1]),:);

n_trial = length(trials_powers);
stimuli_size = zeros(n_trial,n_selected_cells);
for i_trial = 1:n_trial
    for i_cell = 1:n_selected_cells
        distance_center = trials_locations(i_trial,:)-cell_locations(i_cell,:);
        distance_center = round(distance_center + [41 41 91]);
        l23_cells_for_sim(cell_templates(i_cell)).shape(distance_center(1),distance_center(2),distance_center(3));
        stimuli_size(i_trial,i_cell) =trials_powers(i_trial)*...
            l23_cells_for_sim(cell_templates(i_cell)).shape(distance_center(1),distance_center(2),distance_center(3));
    end
end
%%
responses = zeros(n_trial, length(current_template));

responses_cell  = cell([n_selected_cells+1 1]);
for i_cell = 1:(n_selected_cells+1)
    responses_cell{i_cell} = responses;
end
stims = zeros(n_selected_cells,n_trial, length(current_template));

mu_bg = 1/background_rate;
for i_trial = 1:n_trial
    for i_cell = 1:n_selected_cells
        params_sim.V_th= v_th_known(i_cell);
        params_sim.V_reset = v_reset_known(i_cell);
        params_sim.g = [l23_cells_for_sim(cell_templates(i_cell)).g];
        params_sim.gain =  gain_initial(i_cell);
        
        k=stimuli_size(i_trial,i_cell);
        stim = current_template*k;
        stims(i_cell, i_trial,:) = stim;
        [V_vect, spikes]  =lif_glm_sim_v2(stim,params_sim,funcs);
        if delay_params.type==1
        delay_vec=round(normrnd(delay_params.mean,delay_params.std,[sum(spikes) 1]));
        else 
        delay_vec=zeros([sum(spikes) 1]);
        end
        spikes_delay =find(spikes)+delay_vec;
        if spikes_delay < time_max
            if rand(1) < gamma(i_cell)
                responses(i_trial,spikes_delay)=responses(i_trial,spikes_delay)+1;
                responses_cell{i_cell}(i_trial,spikes_delay)=1;
            end
        end
    end
    
    % add background event:
    R = exprnd(mu_bg);
    while R < time_max
        responses(i_trial, max(1,round(R)))=responses(i_trial, max(1,round(R)))+1;
        R = R+exprnd(mu_bg);
    end
end

%% Delay paramters in inferece 
delay_params.mean=35;
delay_params.std=10;
delay_params.type=1;
delay_params.delayed=true;

delay_params.n_grid=200;

% %% Find the minimal phi so that the cell fires at the highest power (direct stim)
% cell_params.gain = max(max(stimuli_size));
% gain_sequence = [1:100]/1e4; %from 0 to 0.01
% n_stimuli_grid_temp=length(gain_sequence);
% 
% cell_params.V_th=v_th_known(1);
% cell_params.V_reset=v_reset_known(i_cell);
% cell_params.g = [l23_cells_for_sim(cell_templates(1)).g];
% cell_params.gain_sd=0.01;
%      
%       
% 
% [Stimuli_grid, Intensity_grid]=Intensity_v8(...
%     gain_sequence,current_template,... % data from exp
%     cell_params,... % estimated parameters
%     funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
%     n_stimuli_grid_temp,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
%     V_threshold,stimulus_threshold,first_only);
% 
% for i_grid = 1:n_stimuli_grid_temp
%     if sum(Intensity_grid{i_grid}) > 0.99
%         gain_lower_bound= gain_sequence(i_grid);
%         break;
%     end
% end
% 
% %% Find the upper bound of gains:
% % The spikes all happen before 20 (~1ms)
% too_early = 20;
% cell_params.gain = max(max(stimuli_size));
% gain_sequence = 0.03+[1:100]/5e2; %from 0.03 to 0.2
% n_stimuli_grid_temp=length(gain_sequence);
% 
% cell_params.V_th=v_th_known(1);
% cell_params.V_reset=v_reset_known(i_cell);
% cell_params.g = [l23_cells_for_sim(cell_templates(1)).g];
% cell_params.gain_sd=0.01;
%      
%      
% [Stimuli_grid, Intensity_grid]=Intensity_v8(...
%     gain_sequence,current_template,... % data from exp
%     cell_params,... % estimated parameters
%     funcs,sd_range, ... % parameters for the GLM model, and range of sd to consider
%     n_stimuli_grid_temp,n_grid_voltage,n_delay_grid,gap_stimuli,... % discrete approximation parameters
%     V_threshold,stimulus_threshold,first_only);
% 
% for i_grid = 1:n_stimuli_grid_temp
%     if sum(Intensity_grid{i_grid}(1:too_early) ) > 0.99
%         gain_upper_bound= gain_sequence(i_grid);
%         break;
%     end
% end
%% Full grid search 
responses_reg=responses;
stims_reg=stims;

loss_type=2;
% Initial values
initial_values = [gain_lower_bound 0.5 gain_lower_bound 0.5];
gain_grid =[0.005 0.008 0.01 0.015 0.02 0.03];
for i = 1:length(gain_grid)
    for j = 1:length(gain_grid)
        initial_values( (i-1)*length(gain_grid)+j,:) = [gain_grid(i) 0.5 gain_grid(j) 0.5];
    end
end
%
in_params.g =  [ l23_cells_for_sim(cell_templates).g];

% LIF-GLM fits
%-------------------------------------%
[stats_conv] = fit_lifglm_v6(responses_reg, stims_reg,in_params,...
    background_rate,v_reset_known,v_th_known,first_only,loss_type,...
    gain_lower_bound,gain_upper_bound, initial_values,delay_params);

% Output:
stats_conv

%% Alternate grid search 
responses_reg=responses;
stims_reg=stims;

loss_type=2;
% Initial values
gain_grid =[0.005 0.008 0.01 0.015 0.02 0.03];
gain_fit = 0.01*ones(n_selected_cells,1);
for i_cell = 1:n_selected_cells
    for i = 1:length(gain_grid)
       initial_values(i,:) = reshape([gain_fit'; 0.5*ones(n_selected_cells,1)'],  [1 2*n_selected_cells]);
       initial_values(i,2*(i_cell-1)+1)=gain_grid(i);
    end
%
in_params.g =  [l23_cells_for_sim(cell_templates).g];

% LIF-GLM fits
%-------------------------------------%
[stats_conv] = fit_lifglm_v6(responses_reg, stims_reg,in_params,...
    background_rate,v_reset_known,v_th_known,first_only,loss_type,...
    gain_lower_bound,initial_values,delay_params);
gain_fit(i_cell)=stats_conv(2*(i_cell-1)+1);
fprintf('%d cell\n',i_cell);
end

initial_values=reshape([gain_fit'; 0.5*ones(n_selected_cells,1)'], [1 2*n_selected_cells]);
[stats_conv] = fit_lifglm_v6(responses_reg, stims_reg,in_params,...
    background_rate,v_reset_known,v_th_known,first_only,loss_type,...
    gain_lower_bound,initial_values,delay_params);
% Output:
stats_conv

%% Feed them the true values 
responses_reg=responses;
stims_reg=stims;
in_params.g =  [l23_cells_for_sim(cell_templates).g];
gain_lower_bound=0.0031;
loss_type=2;

initial_values=reshape([[0.01 0.02]; 1*ones(n_selected_cells,1)'], [1 2*n_selected_cells]);
[stats_conv] = fit_lifglm_v6(responses_reg, stims_reg,in_params,...
    background_rate,v_reset_known,v_th_known,first_only,loss_type,...
    gain_lower_bound,gain_upper_bound,initial_values,delay_params);
% Output:
stats_conv

%% Solve for each cell one by one 
responses_reg=responses;
stims_reg=stims;
in_params.g =  [l23_cells_for_sim(cell_templates).g];
gain_lower_bound=0.0031;
gain_upper_bound = 0.2;
loss_type=2;
single_connection=true;

% Only assigh values in the first column 
gain_grid =[0.01 0.02 0.03];
initial_values = zeros(length(gain_grid),2*n_selected_cells);
initial_values(:,1)=gain_grid;

[stats_conv,~] = fit_lifglm_v6(responses_reg, stims_reg,in_params,...
    background_rate,v_reset_known,v_th_known,first_only,loss_type,...
    gain_lower_bound,gain_upper_bound,initial_values,delay_params,single_connection);
% Output:
stats_conv

%%
xrange=gain_lower_bound + (1:20)/500;
yrange=gain_lower_bound + (1:20)/500;
n_grid=length(current_template);
v_trace=zeros(2, n_trial,n_grid);
for i_cell = 1:2
    for i_trial = 1:n_trial
        v_trace(i_cell,i_trial, 1) = 0;
        for i_t = 2:n_grid
            temp1=reshape(stims(i_cell, i_trial,1:(i_t-1)), [i_t-1,1]);
            temp2=reshape(exp( ((1:(i_t-1)) -i_t+1)*in_params.g(i_cell)), [1,i_t-1]);
            v_trace(i_cell, i_trial, i_t) = temp2*temp1;
        end
    end
end

lklh=zeros(length(xrange),length(yrange));
Xmat=zeros(length(xrange),length(yrange));
Ymat=zeros(length(xrange),length(yrange));

for i = 1:length(xrange)
    for j = 1:length(yrange)
        Xmat(i,j)=xrange(i);
        Ymat(i,j)=yrange(j);
        lklh(i,j)=lif_glm_firstspike_loglikelihood_multi([xrange(i) 1 yrange(j) 1], responses, ...
            v_trace,background_rate,v_th_known,linkfunc,loss_type,delay_params);
    end
    fprintf('%d',i);
end
%%
figure(1)
surf(Xmat,Ymat,lklh)
xlabel('Gain 1')
ylabel('Gain 2')
%%
[~, col_idx] = min(min(lklh));
[~,row_idx] = min(lklh(:,col_idx));
[xrange(row_idx), yrange(col_idx)]
% title('True gain = 0.02')
%xlim([0 0.003])
%     figure(2)
%     plot(xrange,log(lklh(:,2)))
%% Assignments:

soft_true=[];
true_time=[];

soft_bg=[];
bg_time =[];
for i_trial =1:n_trial
    events=find(responses(i_trial,:)>0);
    if length(events)>0
        for i = 1:length(events)
            if true_assignments(i_trial,events(i))==1
                soft_true=[soft_true soft_assignments(i_trial,events(i))];
                true_time=[true_time events(i)];
            else
                soft_bg=[soft_bg soft_assignments(i_trial,events(i))];
                bg_time=[bg_time events(i)];
            end
        end
    end
end

plot(true_time,soft_true,'.','MarkerSize',10)
hold on;
plot(bg_time,soft_bg,'.','MarkerSize',10,'col','r')
hold off;
xlim([0,300])
ylim([0,1])
%%

soft_true=[];
true_time=[];

soft_bg=[];
bg_time =[];
for i_trial =1:n_trial
    events=find(responses_reg(i_trial,:)>0);
    if length(events)>0
        for i = 1:length(events)
            if true_assignments(i_trial,events(i))==1
                soft_true=[soft_true responses_reg(i_trial,events(i))];
                true_time=[true_time events(i)];
            else
                soft_bg=[soft_bg responses_reg(i_trial,events(i))];
                bg_time=[bg_time events(i)];
            end
        end
    end
end

plot(true_time,soft_true,'.','MarkerSize',10)
hold on;
plot(bg_time,soft_bg,'.','MarkerSize',10,'col','r')
hold off;
xlim([0,300])
ylim([0,1])

%% Draw the spike times
figure(2)
plot(sort(spike_time_power{1}),1:length(spike_time_power{1}),...
    '.','col','g','MarkerSize',20)
hold on;
plot(sort(spike_time_power{2}),1:length(spike_time_power{2}),...
    '.','col','b','MarkerSize',20)
hold on;
plot(sort(spike_time_power{3}),1:length(spike_time_power{3}),...
    '.','col','r','MarkerSize',20)
hold off;
ylim([0,num_per_power]);
xlim([0,300]);
line([mean(spike_time_power{1}) mean(spike_time_power{1}) ], [0 20],'col','g')
line([mean(spike_time_power{2}) mean(spike_time_power{2}) ], [0 20],'col','b')
line([mean(spike_time_power{3}) mean(spike_time_power{3}) ], [0 20],'col','r')
xlabel('Time (1/20 ms)');
ylabel('Ranks (irrelevant)');

% saveas(2,'./sorted_time.png');
% saveas(2,'./sorted_time_sim.png');

%%
figure(1)
colors=['g' 'b' 'r'];
t_grid = 1:length(current_template);
for i_outer = 1:3
    avg_spikes{i_outer}=mean(responses( (i_outer-1)*num_per_power +(1:num_per_power),: ));
    %     plot(t_grid,Intensity_grid{i_outer},'col',colors(i_outer),'LineStyle','--','LineWidth',1)
    %     hold on;
    plot(t_grid,avg_spikes{i_outer},'col',colors(i_outer),'LineStyle','--','LineWidth',1)
    hold on;
    avg_spikes_true{i_outer}=mean(responses_raw( (i_outer-1)*num_per_power +(1:num_per_power),: ));
    plot(t_grid,avg_spikes_true{i_outer},'col',colors(i_outer),'LineStyle',':','LineWidth',1)
    hold on;
    
end
xlim([0,300]);
xlabel('Time (1/20 ms)')
ylabel('Intensities')
hold off;
%%

figure(1)
colors=['g' 'b' 'r'];
t_grid = 1:length(current_template);
for i_outer = 1:3
    figure(i_outer)
    %     avg_spikes{i_outer}=mean(responses_reg( (i_outer-1)*num_per_power +(1:num_per_power),: ));
    avg_spikes_true{i_outer}=mean(responses_raw( (i_outer-1)*num_per_power +(1:num_per_power),: ));
    
    plot(t_grid,Intensity_grid{i_outer},'col',colors(i_outer),'LineWidth',1)
    hold on;
    %     plot(t_grid,avg_spikes{i_outer},'col',colors(i_outer),'LineStyle','--','LineWidth',1)
    %     hold on;
    plot(t_grid,avg_spikes_true{i_outer},'col',colors(i_outer),'LineStyle',':','LineWidth',1)
    hold on;
    
    xlim([0,300]);
    xlabel('Time (1/20 ms)')
    ylabel('Intensities')
    hold off;
end
% saveas(3,'./Delay_bg.png');

%%

figure(1)
colors=['r' 'g' 'b'];
t_grid = 1:length(current_template);
for i_outer = 1:3
    plot(t_grid,Intensity_grid{i_outer},'col',colors(i_outer),'LineStyle','--','LineWidth',1)
    hold on;
    
    avg_spikes_true{i_outer}=mean(responses_raw( (i_outer-1)*num_per_power +(1:num_per_power),: ));
    
    plot(t_grid,avg_spikes_true{i_outer},'col',colors(i_outer),'LineStyle',':','LineWidth',1)
    hold on;
    
end
xlim([0,300]);
xlabel('Time (1/20 ms)')
