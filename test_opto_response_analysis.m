
base_cell.start_trial = 1;
base_cell.spike_thresh = 8;
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
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice1_cell1.mat';
base_cell.exclude_trials = [1 2 3];

this_cell = 2; % not great location choices
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice1_cell2.mat';
base_cell.spike_thresh = 10;
base_cell.exclude_trials = [1 4];

this_cell = 3;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice2_cell1.mat';
l23_cells(this_cell).do_vc = 0;
base_cell.exclude_trials = [1 3];

this_cell = 4;
l23_cells(this_cell) = base_cell;% really only looking at first spikes here
l23_cells(this_cell).filename = '4_2_slice2_cell2.mat';
l23_cells(this_cell).spike_thresh = 10;
base_cell.exclude_trials = [1 2 3 4 5];

this_cell = 5; % really only looking at first spikes here
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice2_cell3.mat';
l23_cells(this_cell).spike_thresh = 7;
base_cell.exclude_trials = [];

this_cell = 6; % doublets
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice3_cell1.mat';
base_cell.exclude_trials = [1 4 5];

this_cell = 7; % doublets
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice3_cell2.mat';
l23_cells(this_cell).spike_thresh = 7;
% base_cell.exclude_trials = [1 4 5];

this_cell = 8;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice4_cell1.mat';
l23_cells(this_cell).spike_thresh = 12;

this_cell = 9; % NO SPIKES - CONTROL CELL
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice4_cell2.mat';
l23_cells(this_cell).spike_thresh = 50; % no spikes
l23_cells(this_cell).do_vc = 0;
base_cell.exclude_trials = [1 2 3 4 5];

this_cell = 10;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice1_cell1.mat';
l23_cells(this_cell).start_trial = 3;

this_cell = 11; %doublets
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice1_cell2.mat';
l23_cells(this_cell).start_trial = 3;
l23_cells(this_cell).do_vc = 0;

this_cell = 12; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice1_cell3.mat';
l23_cells(this_cell).start_trial = 2;
l23_cells(this_cell).do_vc = 0;

this_cell = 13; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice2_cell1.mat';
l23_cells(this_cell).start_trial = 3;

this_cell = 14; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice2_cell2.mat';

this_cell = 15; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice3_cell1.mat';
l23_cells(this_cell).spike_thresh = 11;

%% l5 cells

this_cell = 1; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice3_cell2.mat';
l23_cells(this_cell).spike_thresh = 16;

this_cell = 2; 
l5_cells(this_cell) = base_cell;
l5_cells(this_cell).filename = '4_3_slice4_cell3.mat';

this_cell = 3; 
l5_cells(this_cell) = base_cell;
l5_cells(this_cell).filename = '4_3_slice5_cell2.mat';

this_cell = 4; 
l5_cells(this_cell) = base_cell;
l5_cells(this_cell).filename = '4_6_slice2_cell3.mat';
l5_cells(this_cell).spike_thresh = 12;

this_cell = 5; 
l5_cells(this_cell) = base_cell;
l5_cells(this_cell).filename = '4_6_slice3_cell3.mat';
l5_cells(this_cell).hpass_filt = 1;l23_cells(this_cell).type = 'l5pyr';

this_cell = 6; 
l5_cells(this_cell) = base_cell;
l5_cells(this_cell).filename = '4_6_slice3_cell4.mat';
l5_cells(this_cell).hpass_filt = 1;



%%
cell_to_run = [9:15];
% cell_analyzed_bu = cell_analyzed;
% clear cell_analyzed
for ii = 1:length(cell_to_run)
    this_cell = cell_to_run(ii)
    l23_cell_analyzed(this_cell) = analyze_opto_response(l23_cells(this_cell));
    
end

%%
cell_features = [[cell_analyzed.v_th]' [cell_analyzed.gain]' [cell_analyzed.g]'];

figure; [~,ax] = plotmatrix(cell_features);
ylabel(ax(1,1),'v-th')
ylabel(ax(2,1),'optogain')
ylabel(ax(3,1),'g')
title(ax(1,1),'v-th')
title(ax(1,1),'optogain')
title(ax(1,1),'v-th')
title(ax(1,2),'optogain')
title(ax(1,3),'g')


















