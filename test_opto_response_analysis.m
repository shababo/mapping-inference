
base_cell.start_trial = 1;
base_cell.spike_thresh = 5;
base_cell.num_spike_locs = 5;
base_cell.do_cc = 0;
base_cell.do_vc = 1;
base_cell.stim_start = 100;
base_cell.cell_pos = [50 50 0];
base_cell.hpass_filt = 1;
base_cell.exclude_trials = 1;
base_cell.filename = '';
base_cell.type = 'l23pyr';

this_cell = 1;
cell(this_cell) = base_cell;
cell(this_cell).filename = '4_2_slice1_cell1.mat';
base_cell.exclude_trials = [1 2 3];

this_cell = 2; % not great location choices
cell(this_cell) = base_cell;
cell(this_cell).filename = '4_2_slice1_cell2.mat';
base_cell.exclude_trials = [1 4];

this_cell = 3;
cell(this_cell) = base_cell;
cell(this_cell).filename = '4_2_slice2_cell1.mat';
cell(this_cell).do_vc = 0;
base_cell.exclude_trials = [1 3];

this_cell = 4;
cell(this_cell) = base_cell;% really only looking at first spikes here
cell(this_cell).filename = '4_2_slice2_cell2.mat';
cell(this_cell).spike_thresh = 10;
base_cell.exclude_trials = [1 2 3 4 5];

this_cell = 5; % really only looking at first spikes here
cell(this_cell) = base_cell;
cell(this_cell).filename = '4_2_slice2_cell3.mat';
cell(this_cell).spike_thresh = 9;
base_cell.exclude_trials = [];

this_cell = 6; % doublets
cell(this_cell) = base_cell;
cell(this_cell).filename = '4_2_slice3_cell1.mat';
base_cell.exclude_trials = [1 4 5];

this_cell = 7; % doublets
cell(this_cell) = base_cell;
cell(this_cell).filename = '4_2_slice3_cell1.mat';
% base_cell.exclude_trials = [1 4 5];

% 
% this_cell = 15; 
% cell(this_cell) = base_cell;
% cell(this_cell).filename = '4_3_slice3_cell2.mat';
% cell(this_cell).type = 'l5pyr';
% cell(this_cell).spike_thresh = 20;
% 
% this_cell = 16; 
% cell(this_cell) = base_cell;
% cell(this_cell).filename = '4_3_slice4_cell3.mat';
% cell(this_cell).type = 'l5pyr';
% 
% this_cell = 17; 
% cell(this_cell) = base_cell;
% cell(this_cell).filename = '4_3_slice4_cell3.mat';
% cell(this_cell).type = 'l5pyr';
% cell(this_cell).spike_thresh = 12;
% 
% this_cell = 16; 
% cell(this_cell) = base_cell;
% cell(this_cell).filename = '4_6_slice4_cell3.mat';
% cell(this_cell).type = 'l5pyr';
% cell(this_cell).hpass_filt = 1;



%%
cell_to_run = [1:6];
cell_analyzed_bu = cell_analyzed;
clear cell_analyzed
for ii = 1:length(cell_to_run)
    this_cell = cell_to_run(ii)
    cell_analyzed(this_cell) = analyze_opto_response(cell(this_cell));
    
end




















