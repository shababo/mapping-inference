
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
%%
cell_to_run = [1:4];
for i = 1:length(cell_to_run)
    cell_analyzed(cell_to_run(i)) = analyze_opto_response(cell(cell_to_run(i)));
end