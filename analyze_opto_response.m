function this_cell = analyze_opto_response(this_cell)

% this_cell is a struct like this
% this_cell.filename
% this_cell.start_inds... etc

load(this_cell.filename)
% package up and preprocess raw data
start_trial = this_cell.start_trial;
spike_thresh = this_cell.spike_thresh;
ca_num_spike_locs = this_cell.ca_num_spike_locs;
ca_spike_thresh = this_cell.ca_spike_thresh;
do_cc = this_cell.do_cc;
cc_num_spike_locs = this_cell.cc_num_spike_locs;
do_vc = this_cell.do_vc;
stim_start = this_cell.stim_start;
cell_pos = this_cell.cell_pos;
do_hpf = this_cell.hpass_filt;
first_spike_only = this_cell.first_spike;

[this_cell.spike_data, this_cell.voltage_data, this_cell.current_data, this_cell.intrinsics] = ...
    preprocess_opto_response(data,start_trial,spike_thresh,...
    ca_num_spike_locs,ca_spike_thresh,do_cc,cc_num_spike_locs,do_vc,stim_start,cell_pos,first_spike_only);

% do instrinsics
% this_cell.intrinsics = analyze_fi_curve(data,cells.fi_trial); NOT IMPLEMENTED
% at least get V_rest
if isfield(this_cell.intrinsics,'data')
    this_cell.intrinsics.v_rest = median(this_cell.intrinsics.data(1:2*20000));
end
% return
% run lif-glm
downsamp = 1;
if do_cc
    num_spike_locs = cc_num_spike_locs;
    spikes = this_cell.voltage_data;
else
    num_spike_locs = ca_num_spike_locs;
    spikes = this_cell.spike_data;
end
distances = zeros(num_spike_locs,1);

if isfield(this_cell.current_data,'current_shape')
    current_template = this_cell.current_data.current_shape;
else
    load('chrome-template-3ms.mat','template');
    current_template = template(1:1001);
end

count = 1;
for k = 1:num_spike_locs


    spike_times = spikes(k).spike_times;
    powers = spikes(k).powers;
    powers = powers(1:end-1)
    distances(k) = norm(spikes(k).location);
    if do_cc
        this_cell.voltage_data(k).distance = distances(k);
    else
        this_cell.spike_data(k).distance = distances(k);
    end
    for i = 1:length(powers)%-1 % don't do highest power - never useful

        for j = 1:length(spike_times{i})
            responses(count,:) = zeros(1,length(current_template(1:downsamp:end)));
            stims(count,:) = powers(i)*current_template(1:downsamp:end);
            stims_ind(count) = k;
            if ~isempty(spike_times{i}{j})
                % FITTING ONLY FIRST SPIKE!!!
                responses(count,floor(spike_times{i}{j}/downsamp)) = 1;
            end
            count = count + 1;
        end
    end    
end

assignin('base','responses',responses)

% g = [.01 .03 .05 .07 .09 .11 .13]*downsamp;
% g = [.0001 .0005 .001 .005 .01 .05 .1]*downsamp;
% g = [.000001 .000005 .00001 .0005 .0001]*downsamp;
% g = [.001 .002 .003 .004 .005]*downsamp;
% g = [.5 1 5]*downsamp;
g = [.001 .010 .020 .030 .040 .050 .060 .070 .080]*downsamp;
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
num_powers = length(powers);
spike_count_means_glmfit_sim = zeros(num_spike_locs,num_powers);
spike_time_means_glmfit_sim = zeros(num_spike_locs,num_powers);
spike_time_std_glmfit_sim = zeros(num_spike_locs,num_powers);

this_cell.glm_sim.sim_spike_times = cell(num_spike_locs,1);
this_cell.glm_sim.sim_whole_cell = cell(num_spike_locs,1);
% return
for k = 1:num_spike_locs

    k
    this_cell.glm_sim.sim_spike_times{k} = cell(length(powers),1);
    this_cell.glm_sim.sim_whole_cell{k} = cell(length(powers),1);
%     powers = powers;
    params_sim.gain = ...
        this_cell.glm_out(min_ind).glm_result.beta(k+2)*sim_scale;
    for j = 1:length(powers)
        this_cell.glm_sim.sim_spike_times{k}{j} = cell(num_sim_trials,1);
        spike_times = [];
        spike_times_first = [];
        sim_whole_cell = zeros(num_sim_trials,size(stims,2));
        for i = 1:num_sim_trials
            [V_vect, spikes] = lif_glm_sim(stims((j-1)*5+1,:),params_sim,funcs);
            sim_whole_cell(i,:) = V_vect;
            spike_times = [spike_times find(spikes)];
            spike_times_first = [spike_times_first find(spikes,1,'first')];
            this_cell.glm_sim.sim_spike_times{k}{j}{i} = find(spikes);
        end
        this_cell.glm_sim.sim_whole_cell{k}{j} = sim_whole_cell;
        spike_count_means_glmfit_sim(k,j) = length(spike_times)/num_sim_trials;
        spike_time_means_glmfit_sim(k,j) = mean(spike_times_first);
        spike_time_std_glmfit_sim(k,j) = std(spike_times_first);
        
    end    
end
this_cell.glm_sim.spike_count_means = spike_count_means_glmfit_sim;
this_cell.glm_sim.spike_time_means = spike_time_means_glmfit_sim;
this_cell.glm_sim.spike_time_std = spike_time_std_glmfit_sim;
