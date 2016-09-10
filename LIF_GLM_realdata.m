% load 'work_spiking_data_0720.mat'
% load 'current_template.mat'

load('spike_detect_revisit_about_to_fit_0810.mat','all_detection_grids')
load('data/current_template.mat')

%%
trial_length = 1500;
downres_rate = 20;
trial_length_downres = trial_length/downres_rate;

stims_t_norm = -1.0*norm_average_current(1:trial_length);
stims_t_norm_downres = downsample(stims_t_norm,downres_rate);
num_spatial_pos = 121;
powers = [25 50 100];
num_powers = length(powers);
num_repeats = 5;
num_trials = num_spatial_pos*num_repeats*num_powers;


stims_x = [repmat((1:num_spatial_pos)',num_repeats*num_powers,1) [25*ones(num_spatial_pos*num_repeats,1);...
            50*ones(num_spatial_pos*5,1); 100*ones(num_spatial_pos*num_repeats,1)]];




stims_x_value = repmat((1:num_spatial_pos)',num_repeats*num_powers,1);
stims_x_vec = zeros(num_trials,num_spatial_pos);


stims_t_downres = zeros(num_trials,trial_length_downres);

count = 1;

for i = 1:num_powers
    
    for l = 1:num_repeats
        for j = 1:sqrt(num_spatial_pos)
            for k = 1:sqrt(num_spatial_pos)
                stims_t_downres(count,:) = stims_t_norm_downres * powers(i);
                stims_x_vec(count,stims_x_value(count)) = 1;
                count = count + 1;
                
            end
        end
    end
end
%%
num_cells = length(all_detection_grids);

spikes_downres = cell(num_cells,1);
spatial_inits = zeros(11,11,num_cells);
all_detection_grids_downres = cell(size(all_detection_grids));

spike_count_locations = zeros(num_cells,sqrt(num_spatial_pos),sqrt(num_spatial_pos));

for m = 1:num_cells
    

    % stim_t = zeros(size(stims_x,1),1500);
    % % stim_t(100:199) = 1;
    % % stims_t = repmat(stim_t,size(stims_x,1),1);

    spikes_downres{m} = zeros(num_trials,trial_length_downres);
    all_detection_grids_downres{m} = cell(1,3);
    count = 1;
    loc_count = 1;
    for i = 1:3
        all_detection_grids_downres{m}{i} = cell(11,11);
        for l = 1:5
            for j = 1:11
                for k = 1:11
                    if l == 1
                        all_detection_grids_downres{m}{i}{j,k} = zeros(5,75);
                    end
                    if l <= size(all_detection_grids{m}{i}{j,k},1)
                        spike_inds = ceil(find(all_detection_grids{m}{i}{j,k}(l,1:1500))/20);
                        spikes_downres{m}(count,spike_inds) = 1;
                        all_detection_grids_downres{m}{i}{j,k}(l,spike_inds) = 70;
                        if ~isempty(spike_inds)
                            spatial_inits(j,k,m) = spatial_inits(j,k,m) + 1;
                        end
                    else
                        spike_inds = [];
                    end
                    
                    offset = 0;
                    spike_count_locations(m,j-offset,k-offset) = spike_count_locations(m,j-offset,k-offset)...
                        + sum(spikes_downres{m}(count,:));
                    count = count + 1;
                    
                end 
            end
        end
    end

end


%%
% 
cell_choice = 1;
spikes = spikes_downres{cell_choice}';
% spikes = [spikes; zeros(25,num_trials)];

count = 1;
spike_locs = [];
for j = 1:sqrt(num_spatial_pos)
    for k = 1:sqrt(num_spatial_pos)
        
        if spike_count_locations(cell_choice,j,k) > 0
            spike_locs = [spike_locs count];
        end
        count = count + 1;
    end
end



num_trials_good = num_trials - num_repeats*num_powers*(num_spatial_pos - length(spike_locs));
bad_trials = [];
for i = 1:num_trials
    if ~any(stims_x(i,1) == spike_locs)
        bad_trials = [bad_trials i];
    end
end
stims_t_downres_trunc = stims_t_downres;
spikes_trunc = spikes;
stims_t_downres_trunc(bad_trials,:) = [];
spikes_trunc(:,bad_trials) = [];
stims_x_trunc = stims_x;
stims_x_trunc(bad_trials,:) = [];
%%

g_tmp = 1/10;
[expg_hyperpol,expg_rheo,expg_stim]=gconv_multidim(stims_t_downres_trunc',spikes_trunc,g_tmp);

full_stim_mat = zeros(trial_length_downres*num_trials_good,num_spatial_pos);

for i = 1:num_trials_good
    full_stim_mat(1+(i-1)*trial_length_downres:i*trial_length_downres,stims_x_trunc(i,1))...
        = expg_stim(:,i); 
end

stim_scale = 1/100;
full_stim_mat = full_stim_mat*stim_scale;


full_stim_mat_nozero = full_stim_mat(:,spike_locs);




%%

link = @(mu) log(exp(mu)-1);  %link = @(mu) mu + log(1-exp(-mu));
derlink = @(mu) exp(mu)./(exp(mu)-1);
invlink = @(resp) log(1 + exp(resp));
F = {link, derlink, invlink};

[betahat_conv,~,stats_conv]=glmfit([expg_hyperpol(:) full_stim_mat_nozero],spikes_trunc(:),'poisson','link',F);

betahat_conv(1:2)

spatial_filt_fit = zeros(num_spatial_pos,1);
count = 1;
for i = 1:num_spatial_pos
    if any(i == spike_locs)
        spatial_filt_fit(i) = betahat_conv(count+2);
        count = count + 1;
    end
end
figure
imagesc(reshape(spatial_filt_fit,sqrt(num_spatial_pos),sqrt(num_spatial_pos))')
title(['Cell ' num2str(cell_choice) ': Spatial Filter, V_{th} = ' num2str(betahat_conv(1)) ', V_{res} = ' num2str(betahat_conv(2))])

%% output spikes


%%%DEFINE PARAMETERS
dt=1; %time step ms
t_end=75; %total run time ms
V_spike=70; %value to draw a spike to, when cell spikes
tau=10; %membrane time constant [ms]
R_m=10; %membrane resistance [MOhm]
g = 1./tau;

E_L = 0;
V_reset = betahat_conv(2);
V_th = -betahat_conv(1);
spatial_filt = zeros(num_spatial_pos,1);
count = 1;
for i = 1:num_spatial_pos
    if any(i == spike_locs)
        spatial_filt(i) = betahat_conv(count+2);
        count = count + 1;
    end
end

%%%DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS
t_vect=0:dt:t_end; 
V_vect=zeros(1,length(t_vect));
V_plot_vect=zeros(1,length(t_vect));
V_plot_vect2=zeros(1,length(t_vect));


%%
%INTEGRATE THE EQUATION tau*dV/dt = -V + E_L + I_e*R_m
spikes_sim = zeros(t_end,num_trials);
voltages_grid = cell(11,11);
spikes_grid = cell(11,11);
clear I_e_vect_mat

for j = 1:num_trials %loop over different I_Stim values
    
    i=1; %index denoting which element of V is being assigned
        
    V_vect(i)=E_L; %first element of V, i.e. value of V at t=0
    V_plot_vect(i)=V_vect(i); %if no spike, then just plot the actual voltage V
    V_plot_vect2(i)=V_vect(i);
    
    NumSpikes=0; %holds number of spikes that have occurred
    
    tao=exprnd(1);
    lambda(i)=log(exp(V_vect(i)-V_th)+1);
    last_spike = 1;
    
    for t=dt:dt:t_end %loop through values of t in steps of df ms        
        %V_inf = E_L + I_e_vect(i)*R_m;
        %V_vect(i+1) = V_inf + (V_vect(i)-V_inf)*exp(-dt/tau);
        
        V_vect(i+1) = V_vect(i) + ...
            (E_L-V_vect(i) + stims_t_downres(j,i)*spatial_filt(stims_x(j,1)))*g; %Euler's method
        lambda(i+1)=log(exp(V_vect(i+1)-V_th)+1);
        
        %if statement below says what to do if voltage crosses threshold
        if sum(lambda(last_spike+1:i+1))>tao %cell spiked
            V_plot_vect2(i+1) = V_vect(i+1);
            V_vect(i+1)=V_reset; %set voltage back to V_reset
            V_plot_vect(i+1)=V_spike; %set vector that will be plotted to show a spike here
            
            NumSpikes=NumSpikes+1; %add 1 to the total spike count
            spikes_sim(i,j)=1;
            if last_spike == i+1
                disp('two spikes in a row!!')
                return
            end
            last_spike = i+1;
            tao=exprnd(1);
%             lambda(1:i)=0;
%             lambda(i+1)=log(exp(V_vect(i+1)-V_th)+1);
        else %voltage didn't cross threshold so cell does not spike
            V_plot_vect(i+1)=V_vect(i+1); %plot actual voltage
            V_plot_vect2(i+1) = V_vect(i+1);
%             lambda(i+1)=log(exp(V_vect(i+1)-V_th)+1);
        end
        i=i+1; %add 1 to index,corresponding to moving forward 1 time step
    end    
    
    [ind1,ind2] = ind2sub([11 11],stims_x(j,1));
    voltages_grid{ind2,ind1} = [voltages_grid{ind2,ind1}; V_plot_vect];
    spikes_grid{ind2,ind1} = [spikes_grid{ind2,ind1}; spikes_sim(:,j)'];
    %MAKE PLOTS
%     figure(2)
%     subplot(sqrt(num_spatial_pos),sqrt(num_spatial_pos),stims_x(j,1))
%     plot(t_vect,V_plot_vect);
%     hold on
%     plot(t_vect,lambda*100)
%     if (PlotNum==1)
%         title('Voltage vs. time');
%     end
%     if (PlotNum==length(I_Stim_vect))
%         xlabel('Time (ms)');
%     end
%     ylabel('Voltage (mV)');
%     ylim([-100 100])
    
    
    if j ~= num_trials
        clear lambda
        V_vect=zeros(1,length(t_vect));
        V_plot_vect=zeros(1,length(t_vect));
    end
end

%%

voltages_power_grids = cell(3,1);
spikes_grids = cell(3,1);
for i = 1:3
    voltages_power_grids{i} = cell(11,11);
    for j = 1:11
        for k = 1:11
            voltages_power_grids{i}{j,k} = voltages_grid{j,k}((i-1)*5+1:i*5,:);
            spikes_grids{i}{j,k} = spikes_grid{j,k}((i-1)*5+1:i*5,:);
        end
    end
end
%%
figure;compare_trace_stack_grid(voltages_power_grids,5,1,[],0,{'raw','detected events'})
title(['Cell ' num2str(cell_choice) ': Output From LIF-GLM fit'])
drawnow
figure;compare_trace_stack_grid(all_detection_grids_downres{cell_choice},5,1,[],0,{'raw','detected events'})
title(['Cell ' num2str(cell_choice) ': Data (detected spikes)']);

%% summarize events

spikes_per_location_data = zeros(size(voltages_power_grids{1}));
spikes_per_location_sim = zeros(size(voltages_power_grids{1}));
first_spike_latency_data = zeros(size(voltages_power_grids{1}));
first_spike_latency_sim = zeros(size(voltages_power_grids{1}));
spikes_per_location_var_data = zeros(size(voltages_power_grids{1}));
spikes_per_location_var_sim = zeros(size(voltages_power_grids{1}));
first_spike_latency_var_data = zeros(size(voltages_power_grids{1}));
first_spike_latency_var_sim = zeros(size(voltages_power_grids{1}));

for m = 1:length(voltages_power_grids)
    for i = 1:size(voltages_power_grids{m},1)
        for j = 1:size(voltages_power_grids{m},2)
            
            spikes_per_location_sim(i,j) = mean(sum(spikes_grids{m}{i,j}),2);
            spikes_per_location_var_sim(i,j) = var(sum(spikes_grids{m}{i,j}),2);
            
            spikes_per_location_data(i,j) = mean(sum(all_detection_grids_downres{cell_choice}{m}{i,j}),2);
            spikes_per_location_var_data(i,j) = var(sum(all_detection_grids_downres{cell_choice}{m}{i,j}),2);
            
            for k = 1:size(

