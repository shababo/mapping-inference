
base_cell.start_trial = 1;
base_cell.spike_thresh = 7.5;
base_cell.num_spike_locs = 5;
base_cell.do_cc = 0;
base_cell.do_vc = 1;
base_cell.stim_start = 100;
base_cell.cell_pos = [50 50 0];
base_cell.filename = '';

this_cell = 1;
cell(this_cell) = base_cell;
cell(this_cell).filename = '4_2_slice1_cell1.mat';

this_cell = 2;
cell(this_cell) = base_cell;
cell(this_cell).filename = '4_2_slice1_cell2.mat';
cell(this_cell).spike_thresh = 12;

this_cell = 3;
cell(this_cell) = base_cell;
cell(this_cell).filename = '4_2_slice2_cell1.mat';
cell(this_cell).do_vc = 0;

this_cell = 4;
cell(this_cell) = base_cell;
cell(this_cell).filename = '4_2_slice2_cell2.mat';
cell(this_cell).spike_thresh = 24;

this_cell = 5;
cell(this_cell) = base_cell;
cell(this_cell).filename = '4_2_slice2_cell3.mat';
cell(this_cell).spike_thresh = 20;

this_cell = 6; % doublets
cell(this_cell) = base_cell;
cell(this_cell).filename = '4_2_slice3_cell1.mat';
%%
cell_to_run = [1:6];
cell_analyzed_bu = cell_analyzed;
clear cell_analyzed
for ii = 1:length(cell_to_run)
    this_cell = cell_to_run(ii)
    cell_analyzed(this_cell) = analyze_opto_response(cell(this_cell));
    
end




















