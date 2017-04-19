function cell = analyze_opto_response(cell)

% cell is a struct like this
% cell.filename
% cell.start_inds... etc

load(cell.filename)
% package up and preprocess raw data
start_trial = cell.start_trial;
spike_thresh = cell.spike_thresh;
num_spike_locs = cell.num_spike_locs;
do_cc = cell.do_cc;
do_vc = cell.do_vc;
stim_start = cell.stim_start;
cell_pos = cell.cell_pos;
do_hpf = cell.hpass_filt

[cell.spike_data, cell.voltage_data, cell.current_data] = ...
    preprocess_opto_response(data,start_trial,spike_thresh,...
    num_spike_locs,do_cc,do_vc,stim_start,cell_pos);

% do instrinsics
% cell.intrinsics = analyze_fi_curve(data,cells.fi_trial); NOT IMPLEMENTED
% at least get V_rest
if isfield(cell,'instrinsics')
    cell.intrinsics.v_rest = median(cell.intrinsics.data(1:2*20000));
end
    
% run lif-glm
downsamp = 1;
distances = zeros(num_spike_locs,1);

if isfield(cell.current_data,'current_shape')
    current_template = cell.current_data.current_shape;
else
    load('chrome-template-3ms.mat','template');
    current_template = template;
end

count = 1;
for k = 1:num_spike_locs

    spike_times = cell.spike_data(k).spike_times;
    powers = cell.spike_data(k).powers;
    distances(k) = norm(cell.spike_data(k).location);
    cell.spike_data(k).distance = distances(k);
    for i = 1:length(spike_times)

        for j = 1:length(spike_times{i})
            responses(count,:) = zeros(1,length(current_template(1:downsamp:end)));
            stims(count,:) = powers(i)*current_template(1:downsamp:end);
            stims_ind(count) = k;
            if ~isempty(spike_times{i}{j})
                responses(count,floor(spike_times{i}{j}/downsamp)) = 1;
            end
            count = count + 1;
        end
    end    
end

g = [.01 .03 .05 .07 .09 .11 .13]*downsamp;
cell.glm_params.g = g;
cell.glm_params.downsamp = downsamp;

devs = zeros(size(g));
for i = 1:length(g)
    g(i)
    params.g = g(i);
    [cell.glm_out(i).glm_result] = ...
        fit_lifglm(responses,stims,stims_ind,params);
    cell.glm_out(i).dev = cell.glm_out(i).glm_result.dev;
end

[min_dev, min_ind] = min([cell.glm_out.dev]);
cell.g = cell.glm_params.g(min_ind);
this_glm_out = cell.glm_out(min_ind).glm_result;
cell.v_th = this_glm_out.beta(1);
cell.v_reset = this_glm_out.beta(2);
cell.gain = this_glm_out.beta(3);
cell.th_gain_ratio = this_glm_out.beta(1)/this_glm_out.beta(3);

sim_scale = 1;
params_sim.V_th = cell.v_th*sim_scale;
params_sim.V_reset = cell.v_reset*sim_scale;
num_sim_trials = 50;
params_sim.g = cell.g;
funcs.invlink = @invlink_test;
num_locs = length(cell.spike_data(k));
num_powers = length(cell.spike_data(1).powers);
spike_count_means_glmfit_sim = zeros(num_locs,num_powers);
spike_time_means_glmfit_sim = zeros(num_locs,num_powers);
spike_time_std_glmfit_sim = zeros(num_locs,num_powers);


for k = 1:num_locs

    k

    powers = cell.spike_data(k).powers;
    params_sim.gain = ...
        cell.glm_out(min_ind).glm_result.beta(k+2)*sim_scale;
    for j = 1:length(powers)
        spike_times = [];
        sim_whole_cell = zeros(num_sim_trials,size(stims,2));
        for i = 1:num_sim_trials
    
            [V_vect, spikes] = lif_glm_sim(stims((j-1)*5+1,:),params_sim,funcs);
            sim_whole_cell(i,:) = V_vect;
            spike_times = [spike_times find(spikes,1,'first')];
        end
        cell.glm_sim(k).sim_whole_cell = sim_whole_cell;
        spike_count_means_glmfit_sim(k,j) = length(spike_times)/num_sim_trials;
        spike_time_means_glmfit_sim(k,j) = mean(spike_times);
        spike_time_std_glmfit_sim(k,j) = std(spike_times);
        
    end    
end
cell.glm_sim(k).spike_count_means = spike_count_means_glmfit_sim;
cell.glm_sim(k).spike_time_means = spike_time_means_glmfit_sim;
cell.glm_sim(k).spike_time_std = spike_time_std_glmfit_sim;
