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
[cell.spike_data, cell.voltage_data, cell.current_data] = ...
    preprocess_opto_response(data,start_trial,spike_thresh,...
    num_spike_locs,do_cc,do_vc,stim_start);

% do instrinsics
cell.intrinsics = analyze_fi_curve(data,cells.fi_trial);

% run lif-glm
cell.lif_glm_out = fit_lifglm();

