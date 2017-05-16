function this_cell = analyze_opto_response(this_cell)

% this_cell is a struct like this
% this_cell.filename
% this_cell.start_inds... etc

load(this_cell.filename)
% package up and preprocess raw data
start_trial = this_cell.start_trial;
cc_spike_thresh = this_cell.cc_spike_thresh;
ca_num_spike_locs = this_cell.ca_num_spike_locs;
ca_spike_thresh = this_cell.ca_spike_thresh;
do_cc = this_cell.do_cc;
cc_num_spike_locs = this_cell.cc_num_spike_locs;
do_vc = this_cell.do_vc;
stim_start = this_cell.stim_start;
cell_pos = this_cell.cell_pos;
do_hpf = this_cell.hpass_filt;
first_spike_only = this_cell.first_spike;

trial_dur = this_cell.trial_dur;

[this_cell.spike_data, this_cell.voltage_data, this_cell.current_data, this_cell.intrinsics] = ...
    preprocess_opto_response(data,start_trial,cc_spike_thresh,...
    ca_num_spike_locs,ca_spike_thresh,do_cc,cc_num_spike_locs,do_vc,...
    stim_start,cell_pos,first_spike_only,trial_dur,do_hpf);

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

    these_spikes = this_cell.voltage_data;
else
    num_spike_locs = ca_num_spike_locs;
    these_spikes = this_cell.spike_data;
end


% if isfield(this_cell.current_data,'current_shape')
%     current_template = this_cell.current_data.current_shape;
% else
    disp('using general template')
    load('chrome-template-3ms.mat','template');
%     current_template = template(100:1001);
% end
current_template = template(1:trial_dur);
count = 1;
if this_cell.use_shape
    if this_cell.do_vc
        l23_average_shape_norm = ...
            this_cell.current_data.upres_shape;
        upres_x_locs = this_cell.current_data.upres_x_locs;
        upres_y_locs = this_cell.current_data.upres_y_locs;
        upres_z_locs = this_cell.current_data.upres_z_locs;
        load('l23_centered_upres_cell.mat')
        upres_x_locs = -40:1:40;
        upres_y_locs = -40:1:40;
        upres_z_locs = -90:1:90;
%             this_cell.current_data.shape_max/max(this_cell.current_data.shape_max(:));
    else
%         load('l23_mean_shape.mat')
%         load('l23_average_shape_smooth.mat')
%         load('l23_average_shape_centered.mat')
        load('l23_centered_upres_cell.mat')
        upres_x_locs = -40:1:40;
        upres_y_locs = -40:1:40;
        upres_z_locs = -90:1:90;
%         load('l23_gauss_cell_shape.mat')
%         l23_average_shape_norm = l23_gauss_cell_shape;
    end
end
% l23_average_shape_norm = l23_average_shape_norm + .2;
if this_cell.use_shape && ~isempty(this_cell.fit_locs)
    fit_locs = this_cell.fit_locs;
else
    fit_locs = 1:num_spike_locs;
end
num_fit_locs = length(fit_locs);
distances = zeros(num_spike_locs,1);

for kk = 1:num_fit_locs

    k = fit_locs(kk)
    data_spike_times = these_spikes(k).spike_times;
    powers = these_spikes(k).powers;
    powers = powers(1:end-1);%-1 % don't do highest power - never useful
    distances(k) = norm(these_spikes(k).location);

    if do_cc
        this_cell.voltage_data(k).distance = distances(k);
    else
        this_cell.spike_data(k).distance = distances(k);
    end

    
    if this_cell.use_shape
        if k == 1
            these_spikes(k).location = [0 0 0];
            distances(k) = 0;
        end
        this_loc = these_spikes(k).location;
        
        location_ind = [find(this_loc(1) == upres_x_locs) ...
                        find(this_loc(2) == upres_y_locs) ...
                        find(this_loc(3) == upres_z_locs)];
        location_gain = l23_average_shape_norm(location_ind(1),location_ind(2),location_ind(3));
    else
        location_gain = 1;
    end
    for i = 1:length(powers)


        for j = 1:length(data_spike_times{i})
            responses(count,:) = zeros(1,length(current_template(1:downsamp:end)));

            stims(count,:) = powers(i)*current_template(1:downsamp:end)*location_gain;

            if this_cell.use_shape
                stims_ind(count) = 1;
            else
                stims_ind(count) = kk;
            end
            if ~isempty(data_spike_times{i}{j})
                responses(count,floor(data_spike_times{i}{j}/downsamp)) = 1;

            end
            count = count + 1;
        end
    end    
end

assignin('base','responses',responses)
assignin('base','stims',stims)

% g = [.01 .03 .05 .07 .09 .11 .13]*downsamp;
% g = [.0001 .0005 .001 .005 .01 .05 .1]*downsamp;
% g = [.000001 .000005 .00001 .0005 .0001]*downsamp;

g = [.001 .002 .003 .004 .005]*downsamp;
% g = [.025]*downsamp;
% g = [.005 .010 .015 .020 .025 .03]*downsamp; % MAIN ONE

this_cell.glm_params.g = g;
this_cell.glm_params.downsamp = downsamp;

devs = zeros(size(g));
for i = 1:length(g)
%     k
    g(i)
    params.g = g(i);
    [this_cell.glm_out(i).glm_result] = ...
        fit_lifglm(responses,stims,stims_ind,params);
    size(this_cell.glm_out(i).glm_result.beta)
    if this_cell.use_shape
        this_cell.glm_out(i).glm_result.beta(4:num_spike_locs+2) = this_cell.glm_out(i).glm_result.beta(3);
    end
    this_cell.glm_out(i).dev = this_cell.glm_out(i).glm_result.dev;
end

[min_dev, min_ind] = min([this_cell.glm_out.dev]);
this_cell.g = this_cell.glm_params.g(min_ind);
this_glm_out = this_cell.glm_out(min_ind).glm_result;

this_cell.v_th = this_glm_out.beta(1);
this_cell.v_reset = this_glm_out.beta(2);
this_cell.gain = this_glm_out.beta(3:end);
this_cell.th_gain_ratio = this_cell.v_th/this_cell.gain(1);


sim_scale = this_cell.glm_sim_scale;

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
this_cell.glm_sim.lambda = cell(num_spike_locs,1);
% return
lambda = zeros(size(responses));
trial_dur = size(responses,2);
count = 1;
for k = 1:num_spike_locs

    k
    data_spike_times = these_spikes(k).spike_times;
    this_cell.glm_sim.sim_spike_times{k} = cell(length(powers),1);
    this_cell.glm_sim.sim_whole_cell{k} = cell(length(powers),1);
%     powers = powers;
%     location_ind = round([(these_spikes(k).location(1:2)/10 + 5) (these_spikes(k).location(3)/20 + 4)]);
%     location_gain = l23_average_shape_norm(location_ind(1),location_ind(2),location_ind(3));
    params_sim.gain = ...
...%         this_cell.glm_out(min_ind).glm_result.beta(k+2)*sim_scale;
        this_cell.gain(k)*sim_scale;
    if this_cell.use_shape
        
        this_loc = these_spikes(k).location;
        
        location_ind = [find(this_loc(1) == upres_x_locs) ...
                        find(this_loc(2) == upres_y_locs) ...
                        find(this_loc(3) == upres_z_locs)];
                    
        location_gain = l23_average_shape_norm(location_ind(1),location_ind(2),location_ind(3));
    else
        location_gain = 1;
    end
%     these_stims = stims((k-1)*5*num_powers+(1:5*num_powers),:);

    for j = 1:length(powers)
        this_cell.glm_sim.sim_spike_times{k}{j} = cell(num_sim_trials,1);
        sim_spike_times = [];
        spike_times_first = [];
        this_stim  = powers(j)*current_template(1:downsamp:end)*location_gain;
        sim_whole_cell = zeros(num_sim_trials,size(this_stim,2));
        for i = 1:length(data_spike_times{j})
            this_response = responses(count,:);
            [~, ~, lambda(count,:)] = lif_glm_sim(this_stim,params_sim,funcs,this_response);
%             lambda((count-1)*trial_dur+1:count*trial_dur) = this_lambda;
        end
        for i = 1:num_sim_trials
            [V_vect, spikes] = lif_glm_sim(this_stim,params_sim,funcs);
            sim_whole_cell(i,:) = V_vect;
            sim_spike_times = [sim_spike_times find(spikes)];
            spike_times_first = [spike_times_first find(spikes,1,'first')];
            this_cell.glm_sim.sim_spike_times{k}{j}{i} = find(spikes);
        end
        this_cell.glm_sim.sim_whole_cell{k}{j} = sim_whole_cell;
        spike_count_means_glmfit_sim(k,j) = length(sim_spike_times)/num_sim_trials;
        spike_time_means_glmfit_sim(k,j) = mean(spike_times_first);
        spike_time_std_glmfit_sim(k,j) = std(spike_times_first);
        
    end    
end
this_cell.glm_sim.spike_count_means = spike_count_means_glmfit_sim;
this_cell.glm_sim.spike_time_means = spike_time_means_glmfit_sim;
this_cell.glm_sim.spike_time_std = spike_time_std_glmfit_sim;
this_cell.lambda = lambda;

return

if this_cell.do_vc
    
    num_sim_trials = 10;
    test_locs = [-30:10:30; -30:10:30; -60:20:60];
    test_shape_size = size(test_locs);
    this_shape = this_cell.current_data.upres_shape;
    shape_size = size(test_locs);
    
    spike_count_means_glmfit_sim = zeros([test_shape_size num_powers]);
    spike_time_means_glmfit_sim = zeros([test_shape_size num_powers]);
    spike_time_std_glmfit_sim = zeros([test_shape_size num_powers]);

    this_cell.glm_sim.shape_sim.sim_spike_times = cell(test_shape_size);
    this_cell.glm_sim.shape_sim.sim_whole_cell = cell(test_shape_size);

    for k = 1:numel(test_locs)

        k
        this_cell.glm_sim.shape_sim.sim_spike_times{k} = cell(length(powers),1);
        this_cell.glm_sim.shape_sim.sim_whole_cell{k} = cell(length(powers),1);
    %     powers = powers;
%         location_ind = round([(these_spikes(k).location(1:2)/10 + 5) (these_spikes(k).location(3)/20 + 4)]);
        [test_indx, test_indy, test_indz] = lin2sub(shape_size,k);
        this_loc = test_locs(k,:);
        
        location_ind = [find(this_loc(1) == upres_x_locs) ...
                        find(this_loc(2) == upres_y_locs) ...
                        find(this_loc(3) == upres_z_locs)];
                    
        location_gain = l23_average_shape_norm(location_ind(1),location_ind(2),location_ind(3));
        params_sim.gain = ...
    ...%         this_cell.glm_out(min_ind).glm_result.beta(k+2)*sim_scale;
            this_cell.gain(1)*sim_scale;
%         these_stims = stims(1:5*num_powers,:);

        for j = 1:length(powers)
            this_stim  = powers(j)*current_template(1:downsamp:end)*location_gain;
            this_cell.glm_sim.shape_sim.sim_spike_times{k}{j} = cell(num_sim_trials,1);
            data_spike_times = [];
            spike_times_first = [];
            sim_whole_cell = zeros(num_sim_trials,size(these_stims,2));
            for i = 1:num_sim_trials
                [V_vect, spikes] = lif_glm_sim(these_stims((j-1)*5+1,:),params_sim,funcs);
                sim_whole_cell(i,:) = V_vect;
                data_spike_times = [data_spike_times find(spikes)];
                spike_times_first = [spike_times_first find(spikes,1,'first')];
                this_cell.glm_sim.sim_spike_times{k}{j}{i} = find(spikes);
            end
            this_cell.glm_sim.sim_whole_cell{k}{j} = sim_whole_cell;
            spike_count_means_glmfit_sim(k,j) = length(data_spike_times)/num_sim_trials;
            spike_time_means_glmfit_sim(k,j) = mean(spike_times_first);
            spike_time_std_glmfit_sim(k,j) = std(spike_times_first);
        end
    end    
end

this_cell.glm_sim.shape_sim.spike_count_means = spike_count_means_glmfit_sim;
this_cell.glm_sim.shape_sim.spike_time_means = spike_time_means_glmfit_sim;
this_cell.glm_sim.shape_sim.spike_time_std = spike_time_std_glmfit_sim;














