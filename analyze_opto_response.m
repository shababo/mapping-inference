function this_cell = analyze_opto_response(this_cell)

% this_cell is a struct like this
% this_cell.filename
% this_cell.start_inds... etc

load(this_cell.filename)
% package up and preprocess raw data
start_trial = this_cell.start_trial;
spike_thresh = this_cell.spike_thresh;
num_spike_locs = this_cell.num_spike_locs;
do_cc = this_cell.do_cc;
do_vc = this_cell.do_vc;
stim_start = this_cell.stim_start;
cell_pos = this_cell.cell_pos;
do_hpf = this_cell.hpass_filt

[this_cell.spike_data, this_cell.voltage_data, this_cell.current_data, this_cell.intrinsics] = ...
    preprocess_opto_response(data,start_trial,spike_thresh,...
    num_spike_locs,do_cc,do_vc,stim_start,cell_pos);
return
% do instrinsics
% this_cell.intrinsics = analyze_fi_curve(data,cells.fi_trial); NOT IMPLEMENTED
% at least get V_rest
if isfield(this_cell.intrinsics,'data')
    this_cell.intrinsics.v_rest = median(this_cell.intrinsics.data(1:2*20000));
end
    
% run lif-glm
downsamp = 1;
distances = zeros(num_spike_locs,1);

if isfield(this_cell.current_data,'current_shape')
    current_template = this_cell.current_data.current_shape;
else
    load('chrome-template-3ms.mat','template');
    current_template = template(1:1001);
end

count = 1;
for k = 1:num_spike_locs

    spike_times = this_cell.spike_data(k).spike_times;
    powers = this_cell.spike_data(k).powers;
    distances(k) = norm(this_cell.spike_data(k).location);
    this_cell.spike_data(k).distance = distances(k);
    for i = 1:length(powers)%-1 % don't do highest power - never useful

        for j = 1:length(spike_times{i})
            responses(count,:) = zeros(1,length(current_template(1:downsamp:end)));
            stims(count,:) = powers(i)*current_template(1:downsamp:end);
            stims_ind(count) = k;
            if ~isempty(spike_times{i}{j})
                % FITTING ONLY FIRST SPIKE!!!
                responses(count,floor(spike_times{i}{j}(1)/downsamp)) = 1;
            end
            count = count + 1;
        end
    end    
end

assignin('base','responses',responses)

g = [.01 .03 .05 .07 .09 .11 .13]*downsamp;
this_cell.glm_params.g = g;
this_cell.glm_params.downsamp = downsamp;

devs = zeros(size(g));
for i = 1:length(g)
    g(i)
    params.g = g(i);
    [this_cell.glm_out(i).glm_result] = ...
        fit_lifglm(responses,stims,stims_ind,params);
    this_cell.glm_out(i).dev = this_cell.glm_out(i).glm_result.dev;
end

[min_dev, min_ind] = min([this_cell.glm_out.dev]);
this_cell.g = this_cell.glm_params.g(min_ind);
this_glm_out = this_cell.glm_out(min_ind).glm_result;
this_cell.v_th = this_glm_out.beta(1);
this_cell.v_reset = this_glm_out.beta(2);
this_cell.gain = this_glm_out.beta(3);
this_cell.th_gain_ratio = this_glm_out.beta(1)/this_glm_out.beta(3);

sim_scale = 1;
params_sim.V_th = this_cell.v_th*sim_scale;
params_sim.V_reset = this_cell.v_reset*sim_scale;
num_sim_trials = 50;
params_sim.g = this_cell.g;
funcs.invlink = @invlink_test;
num_locs = length(this_cell.spike_data);
num_powers = length(this_cell.spike_data(1).powers);
spike_count_means_glmfit_sim = zeros(num_locs,num_powers);
spike_time_means_glmfit_sim = zeros(num_locs,num_powers);
spike_time_std_glmfit_sim = zeros(num_locs,num_powers);


for k = 1:num_locs

    k

    powers = this_cell.spike_data(k).powers;
    params_sim.gain = ...
        this_cell.glm_out(min_ind).glm_result.beta(k+2)*sim_scale;
    for j = 1:length(powers)
        spike_times = [];
        sim_whole_cell = zeros(num_sim_trials,size(stims,2));
        for i = 1:num_sim_trials
    
            [V_vect, spikes] = lif_glm_sim(stims((j-1)*5+1,:),params_sim,funcs);
            sim_whole_cell(i,:) = V_vect;
            spike_times = [spike_times find(spikes,1,'first')];
        end
        this_cell.glm_sim(k).sim_whole_cell = sim_whole_cell;
        spike_count_means_glmfit_sim(k,j) = length(spike_times)/num_sim_trials;
        spike_time_means_glmfit_sim(k,j) = mean(spike_times);
        spike_time_std_glmfit_sim(k,j) = std(spike_times);
        
    end    
end
this_cell.glm_sim(k).spike_count_means = spike_count_means_glmfit_sim;
this_cell.glm_sim(k).spike_time_means = spike_time_means_glmfit_sim;
this_cell.glm_sim(k).spike_time_std = spike_time_std_glmfit_sim;
