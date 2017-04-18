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

[cell.spike_data, cell.voltage_data, cell.current_data] = ...
    preprocess_opto_response(data,start_trial,spike_thresh,...
    num_spike_locs,do_cc,do_vc,stim_start,cell_pos);

% do instrinsics
% cell.intrinsics = analyze_fi_curve(data,cells.fi_trial);

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

g = [.007 .03 .09 .15 .21]*downsamp;
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

