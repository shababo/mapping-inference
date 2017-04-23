clear l23_cells
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
base_cell.fluor = NaN;

this_cell = 1;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice1_cell1.mat';
l23_cells(this_cell).exclude_trials = [1 2 3];
l23_cells(this_cell).fluor = 4488-3768;

this_cell = 2; % not great location choices
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice1_cell2.mat';
l23_cells(this_cell).spike_thresh = 10;
l23_cells(this_cell).exclude_trials = [1 4];
l23_cells(this_cell).fluor = 2548-2284;

this_cell = 3;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice2_cell1.mat';
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).exclude_trials = [1 3];
l23_cells(this_cell).fluor = 3256-2504;

this_cell = 4;
l23_cells(this_cell) = base_cell;% really only looking at first spikes here
l23_cells(this_cell).filename = '4_2_slice2_cell2.mat';
l23_cells(this_cell).spike_thresh = 10;
l23_cells(this_cell).exclude_trials = [1 2 3 4 5];
l23_cells(this_cell).fluor = 2712-2332;

this_cell = 5; % really only looking at first spikes here
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice2_cell3.mat';
l23_cells(this_cell).spike_thresh = 7;
l23_cells(this_cell).exclude_trials = [];
l23_cells(this_cell).fluor = 2797-2366;

this_cell = 6; % doublets
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice3_cell1.mat';
l23_cells(this_cell).exclude_trials = [1 4 5];
l23_cells(this_cell).fluor = 2653-2437;

this_cell = 7; % doublets
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice3_cell2.mat';
l23_cells(this_cell).spike_thresh = 7;
% base_cell.exclude_trials = [1 4 5];
l23_cells(this_cell).fluor = 2763-2331;

this_cell = 8;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice4_cell1.mat';
l23_cells(this_cell).spike_thresh = 12;
l23_cells(this_cell).fluor = 2804-2411;

this_cell = 9; % NO SPIKES - CONTROL CELL
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice4_cell2.mat';
l23_cells(this_cell).spike_thresh = 50; % no spikes
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).exclude_trials = [1 2 3 4 5];
l23_cells(this_cell).fluor = 2544-2528;

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

this_cell = 16; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice4_cell1.mat';
l23_cells(this_cell).spike_thresh = 7;
l23_cells(this_cell).do_vc = 0;

this_cell = 17; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice4_cell2next.mat';
l23_cells(this_cell).spike_thresh = 7;

this_cell = 18; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice5_cell1next.mat';

this_cell = 19; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice1_cell1next.mat';
l23_cells(this_cell).spike_thresh = 15;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).start_trial = 3;

this_cell = 20; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice1_cell2.mat';
l23_cells(this_cell).spike_thresh = 15;
l23_cells(this_cell).start_trial = 2;

this_cell = 21; % NOT A LOT OF SPIKES
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice1_cell3.mat';
l23_cells(this_cell).spike_thresh = 8;
l23_cells(this_cell).start_trial = 2;

this_cell = 22; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice2_cell1next.mat';
l23_cells(this_cell).spike_thresh = 8;
l23_cells(this_cell).start_trial = 3;

this_cell = 23; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice2_cell2.mat';
l23_cells(this_cell).spike_thresh = 15;
l23_cells(this_cell).start_trial = 2;

this_cell = 24; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice3_cell1next.mat';
l23_cells(this_cell).do_vc = 0;

this_cell = 25; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice3_cell2next.mat';
l23_cells(this_cell).start_trial = 2;



%%
cell_to_run = [1:25];
% cell_analyzed_bu = cell_analyzed;
% clear cell_analyzed
for ii = 1:length(cell_to_run)
    this_cell = cell_to_run(ii)
    l23_cell_analyzed3(this_cell) = analyze_opto_response(l23_cells(this_cell));
    
end


%% l5 cells
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
base_cell.type = 'l5pyr';

this_cell = 1; 

l5_cells(this_cell) = base_cell;
l5_cells(this_cell).filename = '4_3_slice3_cell2.mat';
l5_cells(this_cell).spike_thresh = 16;

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
l5_cells(this_cell).filename = '4_6_slice3_cell3next.mat';
l5_cells(this_cell).hpass_filt = 1;
% l5_cells(this_cell).spike_thresh = 12;
l5_cells(this_cell).do_vc = 0;

this_cell = 6; 
l5_cells(this_cell) = base_cell;
l5_cells(this_cell).filename = '4_6_slice3_cell4next.mat';
l5_cells(this_cell).hpass_filt = 1;
l5_cells(this_cell).spike_thresh = 12;



%%
cell_to_run = [1:6];
% cell_analyzed_bu = cell_analyzed;
% clear cell_analyzed
for ii = 1:length(cell_to_run)
    this_cell = cell_to_run(ii)
    l5_cell_analyzed_shape(this_cell) = analyze_opto_response(l5_cells(this_cell));
    
end

%%
cell_features = [[l23_cell_analyzed2.v_th]' [l23_cell_analyzed2.gain]' [l23_cell_analyzed2.g]'; ...
                 [l5_cell_analyzed2.v_th]' [l5_cell_analyzed2.gain]' [l5_cell_analyzed2.g]'];

groups = ones(size(cell_features,1),1);
groups(end-5:end) = 2;
colors = zeros(length(groups),3);
colors(1:end-6,3) = 1;
colors(end-5:end,2) = 1;

cell_select = 1:length(groups);
% cell_select(union(find(cell_features(:,2) < 0),[12 15 21 24])) = [];

figure; [~,ax] = gplotmatrix(cell_features(cell_select,:),cell_features(cell_select,:),groups(cell_select));
ylabel(ax(1,1),'v-th')
ylabel(ax(2,1),'optogain')
ylabel(ax(3,1),'g')
title(ax(1,1),'v-th')
title(ax(1,1),'optogain')
title(ax(1,1),'v-th')
title(ax(1,2),'optogain')
title(ax(1,3),'g')

ratio = cell_features(cell_select,1)./cell_features(cell_select,2);


figure; scatter(1:length(cell_select),cell_features(cell_select,1)./cell_features(cell_select,2),[],colors(cell_select,:));
ylabel('V-th to optogain ratio')
xlabel('cell id')

figure
scatter(2.0*ones(length(ratio(groups == 1)),1),ratio(groups == 1),'b' ,'jitter','on', 'jitterAmount',0.05);
hold on
scatter(3.0*ones(length(ratio(groups == 2)),1),ratio(groups == 2),'g' ,'jitter','on', 'jitterAmount',0.05);
xlim([1.5 3.5])
ylim([0 1750])
ylabel('optical sensitivity (expressoin x rheobase)')
set(gca,'xtick',[2 3])
set(gca,'xticklabels',{'l23pyr','l5pyr'})

%% shapes
figure
cells_to_plot = find([l5_cells.do_vc]);
z_slices = size(l5_cell_analyzed_shape(1).current_data.shape_max,3);
l5_average_shape = zeros(9,9,7);
z_depths = {'-90','-50','-20','0','20','50','90'};
count = 1;
for i = 1:length(cells_to_plot)
    cell_i = cells_to_plot(i);
    this_shape = l5_cell_analyzed_shape(cell_i).current_data.shape_max;
    this_shape = this_shape/max(this_shape(:));
    l5_average_shape = l5_average_shape + this_shape;
    for j = 1:z_slices
        subplot(z_slices,length(cells_to_plot)+1,(j-1)*(length(cells_to_plot)+1) + i + 1)
        imagesc(this_shape(:,:,j))
        caxis([0 1])
        axis off
        count = count + 1;
        if j == 1
            title(['cell ' num2str(i)])
        end
    end
end

l5_average_shape = l5_average_shape/length(cells_to_plot);
for j = 1:z_slices
    subplot(z_slices,length(cells_to_plot)+1,(j-1)*(length(cells_to_plot)+1) + 1)
    imagesc(l23_average_shape(:,:,j))
%     axis off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    caxis([0 1])
    
    if j == 1
        title('mean cell')
    end
    ylabel(['z = ' z_depths{j}])
end

%%

figure
cells_to_plot = find([l23_cells.do_vc]);
avg_template = [];
for i = 1:length(cells_to_plot)
    ii = cells_to_plot(i);
    if i == 1
        avg_template = l23_cell_analyzed3(ii).current_data.current_shape;
    else
        avg_template = avg_template + l23_cell_analyzed3(ii).current_data.current_shape;
    end
    plot((1:length(l23_cell_analyzed3(ii).current_data.current_shape))/20000,...
        l23_cell_analyzed3(ii).current_data.current_shape,'color',[.6 .6 .6])
    hold on
end
avg_template = avg_template/length(cells_to_plot);
plot((1:length(avg_template))/20000,avg_template,'k','linewidth',3)
hold off
xlabel('time (sec)')
ylabel('current (a.u.)')
title('est. chrome current shape - 3 msec pulse - l23pyr')
xlim([0 .05])



        















