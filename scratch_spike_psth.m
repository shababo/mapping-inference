%%

load('work_spiking_data_0719.mat')

%% spike detection on grid
 
trace_grids_3_31_s2c2_r2_3 = {traces_by_location_3_31_s2c2_r3_5mw, traces_by_location_3_31_s2c2_r3_10mw, traces_by_location_3_31_s2c2_r3_15mw,...
    traces_by_location_3_31_s2c2_r3_25mw, traces_by_location_3_31_s2c2_r3_50mw, traces_by_location_3_31_s2c2_r3_100mw};

trace_grids_3_29_s1c2_r2 = {traces_by_location_3_29_s1c2_r2_25mw, traces_by_location_3_29_s1c2_r2_50mw, traces_by_location_3_29_s1c2_r2_100mw};

trace_grids_3_29_s1c4_r2 = {traces_by_location_3_29_s1c4_r2_25mw, traces_by_location_3_29_s1c4_r2_50mw, traces_by_location_3_29_s1c4_r2_100mw};

%trace_grids_3_31_s1c1_r4_5 = {traces_by_location_3_31_s1c1_r4_5_25mw, traces_by_location_3_31_s1c1_r4_5_50mw, traces_by_location_3_31_s1c1_r4_%5_100mw};

trace_grids_3_31_s1c2_r4_5 = {traces_by_location_3_31_s1c2_r5_5mw, traces_by_location_3_31_s1c2_r5_10mw, traces_by_location_3_31_s1c2_r5_15mw,...
    traces_by_location_3_31_s1c2_r4_25mw, traces_by_location_3_31_s1c2_r4_50mw, traces_by_location_3_31_s1c2_r4_100mw};

trace_grids_4_5_s2c1_r5 = {traces_by_location_4_5_s2c1_r5_25mw, traces_by_location_4_5_s2c1_r5_50mw, traces_by_location_4_5_s2c1_r5_100mw};

trace_grids_4_6_s3c2_r1 = {traces_by_location_4_6_s3c2_r1_25mw, traces_by_location_4_6_s3c2_r1_50mw, traces_by_location_4_6_s3c2_r1_100mw};

trace_grids_4_6_s3c5_r1 = {traces_by_location_4_6_s3c5_r1_25mw, traces_by_location_4_6_s3c5_r1_50mw, traces_by_location_4_6_s3c5_r1_100mw};

trace_grids_4_6_s3c7_r2 = {traces_by_location_4_6_s3c7_r2_25mw, traces_by_location_4_6_s3c7_r2_50mw, traces_by_location_4_6_s3c7_r2_100mw};

trace_grids_4_6_s3c8_r3 = {traces_by_location_4_6_s3c8_r3_25mw, traces_by_location_4_6_s3c8_r3_50mw, traces_by_location_4_6_s3c8_r3_100mw};

%%

trace_grids = trace_grids_4_6_s3c2_r1;

detection_grids = cell(size(trace_grids));

for i = 1:length(trace_grids)
    
    trace_grid_tmp = trace_grids{i};
    [traces_tmp, rebuild_map] = stack_traces(trace_grid_tmp);

    detection_results = detect_peaks(-1.0*bsxfun(@minus,traces_tmp,median(traces_tmp(:,1300:1500),2)),1.0,20,1,1,0)*70;
    detection_grids{i} = unstack_traces(detection_results,rebuild_map);
    
end
    
detection_results_4_6_s3c2_r1 = detection_results;
detection_grids_4_6_s3c2_r1 = detection_grids;


figure; compare_trace_stack_grid({trace_grids{:},detection_grids_4_6_s3c2_r1{:}},...
    5,1,[],0,{'25 mW', '50 mW', '100 mW'},2)

%%

all_detection_grids = {detection_grids_3_31_s2c2_r2_3(4:end), detection_grids_3_29_s1c2_r2, ...
    detection_grids_3_29_s1c4_r2, detection_grids_3_31_s1c1_r4_5, detection_grids_3_31_s1c2_r4_5(4:end), ...
    detection_grids_4_5_s2c1_r5, detection_grids_4_6_s3c2_r1, detection_grids_4_6_s3c5_r1, ...
    detection_grids_4_6_s3c7_r2, detection_grids_4_6_s3c8_r3};

%%

figure; compare_trace_stack_grid({trace_grids_3_31_s1c1_r4_5{:},detection_grids_3_31_s1c1_r4_5{:}},...
    5,1,[],0,{'25 mW', '50 mW', '100 mW'},2)

%%

psth_grids = cell(3,1);
single_cell_psth_grids = cell(length(all_detection_grids),3);
for i = 1:3
    i
    psth_grids{i} = cell(11,11);
    for  m = 1:length(all_detection_grids)
        m
        single_cell_psth_grids{m,i} = cell(11,11);
        for j = 1:11
            j
            for k = 1:11
                k
                single_cell_psth_grids{m,i}{j,k} =  smoothts(sum(all_detection_grids{m}{i}{j,k}(:,1:1500),1)/size(all_detection_grids{m}{i}{j,k},1),'g',500,100)/2;
                if m == 1
                   psth_grids{i}{j,k} = sum(all_detection_grids{m}{i}{j,k}(:,1:1500),1)/size(all_detection_grids{m}{i}{j,k},1)/length(all_detection_grids);
                else
                   size(psth_grids{i}{j,k})
                   size(sum(all_detection_grids{m}{i}{j,k},1)/size(all_detection_grids{m}{i}{j,k},1))
                   psth_grids{i}{j,k} = psth_grids{i}{j,k} + sum(all_detection_grids{m}{i}{j,k}(:,1:1500),1)/size(all_detection_grids{m}{i}{j,k},1)/length(all_detection_grids);
                end
               
                if m == length(all_detection_grids)
                    psth_grids{i}{j,k} = smoothts(psth_grids{i}{j,k},'g',100,20);
                end
            end
        end
    end
end

%%
figure; compare_trace_stack_grid(psth_grids(1),...
    1,1,[],0,{'25 mW', '50 mW', '100 mW'},1)

%%
for i = 1:1%size(single_cell_psth_grids,1)
    figure
    compare_trace_stack_grid(single_cell_psth_grids(i,:),...
    1,1,[],0,{'25 mW', '50 mW', '100 mW'},1)
    set(gcf,'position',[76 484 1844 510])
end
    


%% make peak and delay images
peak_image = zeros(11,11,3);
delay_image = zeros(11,11,3);

figure
for i = 1:3
    for j = 1:11
        for k = 1:11
            [peak_image(j,k,i) delay_image(j,k,i)] = max(psth_grids{i}{j,k});
        end
    end
    subplot(2,3,i)
    imagesc(peak_image(:,:,i))
    caxis([0 2])
    colorbar
    subplot(2,3,i+3)
    imagesc(delay_image(:,:,i))
    caxis([0 1500])
    colorbar
end

colormap hot


%% gen fake version - must run gendata_fullmodel.m first to get variables in ws

%% stim paramters

% covariance of point spread function
A = diag([500, 500, 1000]);

evoked_params.stim_tau_rise = .0015*20000; % values for chr2 from lin et al 2009 (biophysics)
evoked_params.stim_tau_fall = .013*20000;
evoked_params.stim_amp = 0;
evoked_params.stim_start = .005*20000;
evoked_params.stim_duration = .005*20000;

% effect on postsynaptic cell
evoked_params.stim_tau_rise = .0015*20000; % values for chr2 from lin et al 2009 (biophysics)
evoked_params.stim_tau_fall = .013*20000;
evoked_params.stim_amp = 0;

%% select a postsyanptic cell
cell_layer = 5; % 5A
num_cell_layer_neurons = size(neuron_locations{cell_layer},1);
postsyn_position = neuron_locations{cell_layer}(randi(num_cell_layer_neurons),:);

%% Generate some data

num_repeats = 5;
N = 11*11*num_repeats;
R = 1;

% these parameters govern the time delay, as a function of the
% point-spread-function stimuli for a particular trial
% in seconds
d_mean0 = .000;
d_sigma0 = .002;
d_mean_coef = .005;
d_sigma_coef = .050;

% the neurons being stimulated at each trial
% Z = false(N,K);
Z = zeros(N/num_repeats,3);
trial_grid_locations = zeros(N/num_repeats,2);
count = 1;
for i = 1:11
    for j = 1:11
        
        trial_grid_locations(count,:) = [i j];
        count = count + 1;
        Z((i-1)*11 + j,1) = (i-1)*10 - 50 + postsyn_position(1);
        Z((i-1)*11 + j,2) = (j-1)*10 - 50 + postsyn_position(2);
        Z((i-1)*11 + j,3) = postsyn_position(3);
        
    end
end

Z = repmat(Z,num_repeats,1);
trial_grid_locations = repmat(trial_grid_locations,num_repeats,1);

% probability of firing
% pi_nk = zeros(N,K);
K = 10;
pi_kr = exp(-0.5*squareform(pdist([Z; repmat(postsyn_position,K,1)],'mahalanobis',A)).^2);
pi_nk = pi_kr(1:N,N+1:N+K);
pi_nk_spike = pi_nk;
pi_nk_spike(pi_nk_spike > .5) = 1;
% pi_nk = pi_nk/.5;

%%

% firing delay means and variances
d_mean_nk = d_mean0 + (1.5 - pi_nk).^1*d_mean_coef;
d_sigma_nk = d_sigma0 + (1 - pi_nk)*d_sigma_coef;

% sample "ground truth" firing delay
D = exprnd(d_mean_nk/data_params.dt) + evoked_params.stim_start;
% D(D < evoked_params.stim_start + .002) = evoked_params.stim_start + .002;

% sample "ground truth" stimulations

X = rand(N,K) < pi_nk_spike; %.2 
X(D > 1500) = 0;

all_rasters = cell(K,1);

for k = 1:K
    all_rasters{k} = zeros(size(D,1),1500);
    for  i = 1:size(D,1)
        all_rasters{k}(i,ceil(D(i,k))) = X(i,k);
    end
    all_rasters{k} = unstack_traces(all_rasters{k},trial_grid_locations);
end
%%

psth_grids_fake = cell(11,11);
for  m = 1:length(all_rasters)
        m
        
    for j = 1:11
        j
        for k = 1:11
            k
            if m == 1
               psth_grids_fake{j,k} = sum(all_rasters{m}{j,k}(:,1:1500),1)/size(all_rasters{m}{j,k},1)/length(all_rasters);
            else
               size(psth_grids_fake{j,k})
               size(sum(all_rasters{m}{j,k},1)/size(all_rasters{m}{j,k},1))
               psth_grids_fake{j,k} = psth_grids_fake{j,k} + sum(all_rasters{m}{j,k}(:,1:1500),1)/size(all_rasters{m}{j,k},1)/length(all_rasters);
            end

            if m == length(all_rasters)
                psth_grids_fake{j,k} = smoothts(psth_grids_fake{j,k},'g',100,20)*40;
            end
        end
    end
end
        

figure; compare_trace_stack_grid({psth_grids_fake},...
    1,1,[],0,{'25 mW', '50 mW', '100 mW'},1)
