%% Toy example of new simulations with LIF 
% Generate the postsynaptic event given a stimulus
%% Loading functions and Data generation
clear;
addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Parameters 
rng(12242,'twister');

% load parameters for the model
% run('./Parameters_setup_ground_truth.m')
run('./Parameters_setup_LIF.m')

num_dense=40; % Number of grids
Aval = 400;
A = diag([Aval, Aval, 1500]); %<--- 

run('Parameters_setup_experiment.m')
num_sources = 4;  % Number of locations to stimulate in each trial

bg_params.firing_rate = 20; % background firing rate
%% LIF related Parameters 
% 

% Fix the stimulation to be 100 for now
num_I_Stim=1; 
% 1: 100; 6: 50; 11: 25;    

g=0.1; %membrane time constant [ms]
k_offset = 0.04; % The spatial mark; larger value -> more spikes
% NOTE:
% The values for the V_threshold, V_reset, V_resting are set in lines 96 --
% 106 in Parameters_setup_LIF.m
 

%% Conduct one trial:
grid_index = 1:size(pi_dense_all,2);
K_z = size(pi_dense_local,1);

locations_next = randsample(grid_index, num_sources);

X_next_all = min(.95,sum(pi_dense_all(:,locations_next),2));


evoked_k = k_offset.* X_next_all; %<---- 

% Draw presynaptic events:
presynaptic_events = cell(size(evoked_k,1),1);
for i_cell = 1:size(evoked_k,1)
    if all_amplitudes(i_cell) > 0
        V_th = all_V_th(i_cell);
        V_reset = all_V_reset(i_cell);
        E_L=-28; %resting membrane potential [mV]
        
        k = evoked_k(i_cell);
        presynaptic_events{i_cell} = [];
        if k > 0.00001
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
                %if statement below says what to do if voltage crosses threshold
                if (V_vect(i+1)>V_th) %cell spiked
                    V_vect(i+1)=V_reset; %set voltage back to V_reset
                    presynaptic_events{i_cell} = [presynaptic_events{i_cell} ti];
                    %t
                end
                i=i+1; %add 1 to index,corresponding to moving forward 1 time step
            end
        end
    end
end

% Draw spontaneous events:
background_events = [];
mu_bg = 1000/bg_params.firing_rate;
T_remain =data_params.T;
R = exprnd(mu_bg);
while R < T_remain
    background_events=[background_events R];
    T_remain = T_remain - R;
    R = exprnd(mu_bg);
end


% Merging the presynaptic events
mpp_n = struct();
if isempty(background_events) ==  0
    mpp_n.event_times = background_events;
    mpp_n.amplitudes = rand([1 length(background_events)])...
        .*(bg_params.a_max -bg_params.a_min)+bg_params.a_min;
    mpp_n.assignments = zeros(size(background_events));
else
    mpp_n.event_times=[];
    mpp_n.amplitudes=[];
    mpp_n.assignments=[];
    
end

for i_cell = 1:size(evoked_k,1)
    if all_amplitudes(i_cell) > 0
        if isempty(presynaptic_events{i_cell}) ==  0
            for i = 1:length(presynaptic_events{i_cell})
                if rand(1) > evoked_params.failure_prob
                    mpp_n.event_times = [mpp_n.event_times presynaptic_events{i_cell}(i)];
                    mpp_n.amplitudes = [mpp_n.amplitudes normrnd(all_amplitudes(i_cell),evoked_params.sigma_a)];
                    mpp_n.assignments = [mpp_n.assignments i_cell];
                end
            end
        end
    end
end


% Plotting:

% Stimulus and spiked presynaptic cells

figure(1)

for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*35);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,0.3);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})

xlim([20,460]);
ylim([-900,-400]);


stimulus= scatter(Z_dense(locations_next,1), -Z_dense(locations_next,2),...
    50, 'filled','d');
set(stimulus,'MarkerFaceColor','b');
alpha(stimulus,0.8);
hold on

if sum(mpp_n.assignments > 0) >0
    non_bg = unique( mpp_n.assignments(mpp_n.assignments>0));
    spiked = scatter(all_neuron_locations(non_bg,1), -all_neuron_locations(non_bg,2),...
        20, 'filled','o');
    set(spiked,'MarkerFaceColor','g');
    alpha(spiked,0.8);
    hold on
end

hold off


% Draw the spike train:
figure(2)

ind_vec = mpp_n.assignments  == 0;
x = [mpp_n.event_times(ind_vec); mpp_n.event_times(ind_vec)];
y = [zeros(length(mpp_n.amplitudes(ind_vec)),1) mpp_n.amplitudes(ind_vec)'];

plot(x,y','col','k','Linewidth',4)
hold on;

ind_vec = mpp_n.assignments  > 0;
x = [mpp_n.event_times(ind_vec); mpp_n.event_times(ind_vec)];
y = [zeros(length(mpp_n.amplitudes(ind_vec)),1) mpp_n.amplitudes(ind_vec)'];
plot(x,y','col','g','Linewidth',4)
hold on;

xlim([0,75]);
ylim([0,max(mpp_n.amplitudes)]);

xlabel('Time (ms)');
ylabel('Event size');

hold off

%%
m = [11    44     2     9;   11    44     5     8;   2     1     6    11;   2     1    10     3];

x=[m(:,1) m(:,3)];
y=[m(:,2) m(:,4)];
plot(x',y')

