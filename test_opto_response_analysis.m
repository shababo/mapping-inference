clear l23_cells
base_cell.start_trial = 1;
base_cell.spike_thresh = 8;
base_cell.ca_num_spike_locs = 5;
base_cell.ca_spike_thresh = -10;
base_cell.do_cc = 0;
base_cell.cc_num_spike_locs = 0;
base_cell.do_vc = 1;
base_cell.stim_start = 100;
base_cell.cell_pos = [50 50 0];
base_cell.hpass_filt = 1;
base_cell.exclude_trials = 1;
base_cell.filename = '';
base_cell.type = 'l23pyr';
base_cell.fluor = NaN;
base_cell.first_spike = 0;
base_cell.glm_sim_scale = 1;
base_cell.trial_dur = .030*20000;

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
l23_cells(this_cell).fluor = 3004-2788;

this_cell = 20; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice1_cell2.mat';
l23_cells(this_cell).spike_thresh = 15;
l23_cells(this_cell).start_trial = 2;
l23_cells(this_cell).fluor = 292-189;

this_cell = 21; % NOT A LOT OF SPIKES
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice1_cell3.mat';
l23_cells(this_cell).spike_thresh = 8;
l23_cells(this_cell).start_trial = 2;
l23_cells(this_cell).fluor = 261-202;

this_cell = 22; % START FLUOR HERE??
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice2_cell1next.mat';
l23_cells(this_cell).spike_thresh = 8;
l23_cells(this_cell).start_trial = 3;
l23_cells(this_cell).fluor = 590-220;

this_cell = 23; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice2_cell2.mat';
l23_cells(this_cell).spike_thresh = 15;
l23_cells(this_cell).start_trial = 2;
l23_cells(this_cell).fluor = 720-223;

this_cell = 24; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice3_cell1next.mat';
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).fluor = 878-197;

this_cell = 25; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice3_cell2next.mat';
l23_cells(this_cell).start_trial = 2;
l23_cells(this_cell).fluor = 857-217;

this_cell = 26;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_21_slice2_cell2.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 10;
l23_cells(this_cell).do_vc = 1;

this_cell = 27;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '3_30_slice2_cell1.mat';
l23_cells(this_cell).ca_num_spike_locs = 10;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 10;
l23_cells(this_cell).spike_thresh = 20;

this_cell = 28;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '3_30_slice2_cell2.mat';
l23_cells(this_cell).ca_num_spike_locs = 10;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 10;
l23_cells(this_cell).spike_thresh = 20;

this_cell = 29;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '3_30_slice2_cell3.mat';
l23_cells(this_cell).ca_num_spike_locs = 10;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 10;
l23_cells(this_cell).spike_thresh = 20;

this_cell = 30;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '3_30_slice4_cell2.mat';
l23_cells(this_cell).ca_num_spike_locs = 10;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 10;
l23_cells(this_cell).spike_thresh = 20;

this_cell = 31;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_21_slice3_cell2.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 10;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).spike_thresh = 20;
l23_cells(this_cell).ca_spike_thresh = -15;

this_cell = 32;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_25_slice1_cell1.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 6;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).spike_thresh = 20;
l23_cells(this_cell).ca_spike_thresh = -25;

this_cell = 33;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_25_slice1_cell2.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 6;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).spike_thresh = 20;
l23_cells(this_cell).ca_spike_thresh = -25;

this_cell = 34;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_25_slice2_cell1.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 6;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).spike_thresh = 20;

this_cell = 35;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_25_slice2_cell2.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 6;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).spike_thresh = 20;

this_cell = 36;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_25_slice2_cell3.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 6;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).spike_thresh = 20;

this_cell = 37;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_25_slice3_cell1.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 6;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).spike_thresh = 20;

this_cell = 38;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_25_slice3_cell2.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 6;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).spike_thresh = 20;


%%

cell_to_run = [26:38];
% cell_analyzed_bu = cell_analyzed;
% clear l23_cell_analyzed13 <-- multispike, 14 is one spike

for ii = 1:length(cell_to_run)
    this_cell = cell_to_run(ii)
    l23_cell_analyzed_fulldata_longer(this_cell) = analyze_opto_response(l23_cells(this_cell));

    
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
cells_to_plot = find([l23_cells.do_vc]);
cells_to_plot = setdiff(cells_to_plot,[5 6 27:30]);
z_slices = size(l23_cell_analyzed_preprocessonly(1).current_data.shape_max,3);
l23_average_shape = zeros(9,9,7);
z_depths = {'-90','-50','-20','0','20','50','90'};
count = 1;
for i = 1:length(cells_to_plot)
    cell_i = cells_to_plot(i);
    this_shape = l23_cell_analyzed_preprocessonly(cell_i).current_data.shape_max;
    this_shape = this_shape/max(this_shape(:));
    l23_average_shape = l23_average_shape + this_shape;
    for j = 1:z_slices
        subplot(z_slices,length(cells_to_plot)+1,(j-1)*(length(cells_to_plot)+1) + i + 1)
        imagesc(this_shape(:,:,j))
        caxis([0 1])
        axis off
        count = count + 1;
        if j == 1
            title(['cell ' num2str(cell_i)])
        end
    end
end

l23_average_shape = l23_average_shape/max(l23_average_shape(:));
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
cells_to_plot = setdiff(cells_to_plot,5);
avg_template = [];
colors = parula(100);
svals = zeros(length(cells_to_plot),1);
for i = 1:length(cells_to_plot)
    ii = cells_to_plot(i);
    svals(i) = l23_cell_analyzed3(ii).current_data.current_shape_sv;
end
% svals = floor(svals/max(svals)*100);
svals = floor(svals*100);
for i = 1:length(cells_to_plot)
    ii = cells_to_plot(i);
    if i == 1
        avg_template = l23_cell_analyzed3(ii).current_data.current_shape;
    else
        avg_template = avg_template + l23_cell_analyzed3(ii).current_data.current_shape;
    end
    plot((1:length(l23_cell_analyzed3(ii).current_data.current_shape))/20000,...
        l23_cell_analyzed3(ii).current_data.current_shape,'color',colors(svals(i),:))
    hold on
end
avg_template = avg_template/length(cells_to_plot);
plot((1:length(avg_template))/20000,avg_template,'k','linewidth',3)
hold off
xlabel('time (sec)')
ylabel('current (a.u.)')
title('est. chrome current shape - 3 msec pulse - l23pyr')
xlim([0 .05])


%%

% analyzed_cells = l23_cell_analyzed_fulldata_longer;
cell_select = 26:38;

opt_rheobase = zeros(size(cell_select));
mean_time = zeros(size(cell_select));

for i = 1:length(cell_select)
    
    this_cell = cell_select(i);
    op_rheobase(i) = l23_cell_analyzed_fulldata_longer(this_cell).th_gain_ratio;
    mean_time(i) = l23_cell_analyzed6(this_cell).voltage_data(1).spike_times_means(end);
    
end

% figure; scatter(mean_time,op_rheobase)

figure; plot(mean_time)

%% devs

figure;
cell_select = [26:38];
g = [.0001 .0005 .001 .005 .01 .05 .1];
g = [1.000000000000000e-03,0.002000000000000,0.003000000000000,0.004000000000000,0.005000000000000];
colors = parula(length(cell_select));
for i = 1:length(cell_select)
    
    this_cell = cell_select(i);
%     semilogx(g,[l23_cell_analyzed5(this_cell).glm_out.dev l23_cell_analyzed_fulldata_longer(this_cell).glm_out.dev l23_cell_analyzed6(this_cell).glm_out.dev],'color',colors(i,:))
plot(l23_cell_analyzed_fulldata_longer(this_cell).glm_params.g,[l23_cell_analyzed_fulldata_longer(this_cell).glm_out.dev],'color',colors(i,:))
    hold on
    
end


%% num_spikes

figure;
cell_select = 26:38;
cell_select = setdiff(cell_select,[29 31 32]);
colors = parula(length(cell_select));
count = 1;
for i = 1:length(cell_select)
    this_cell = cell_select(i);
%     subplot(ceil(length(cell_select)/4),4,i)
    if l23_cells(this_cell).do_cc
        for j = 1:length(l23_cell_analyzed_fulldata_longer(this_cell).voltage_data)
            subplot(length(cell_select),10,(i-1)*10 + j)
            plot(l23_cell_analyzed_fulldata_longer(this_cell).voltage_data(j).powers(1:size(l23_cell_analyzed_fulldata_longer(this_cell).glm_sim.spike_count_means,2)),...
                l23_cell_analyzed_fulldata_longer(this_cell).voltage_data(j).num_spike_means(1:size(l23_cell_analyzed_fulldata_longer(this_cell).glm_sim.spike_count_means,2)),'k.-')

            hold on
            plot(l23_cell_analyzed_fulldata_longer(this_cell).voltage_data(j).powers(1:size(l23_cell_analyzed_fulldata_longer(this_cell).glm_sim.spike_count_means,2)),...
                l23_cell_analyzed_fulldata_longer(this_cell).glm_sim.spike_count_means(j,:),'bo-')
%             hold on
%              plot(l23_cell_analyzed21(this_cell).voltage_data(j).powers(1:size(l23_cell_analyzed21(this_cell).glm_sim.spike_count_means,2)),...
%                 l23_cell_analyzed21(this_cell).glm_sim.spike_count_means(j,:),'m')
            
            count = count + 1;
            title(mat2str(l23_cell_analyzed_fulldata_longer(this_cell).voltage_data(j).location))
            ylim([0 3.0])
            xlim([0 110])
            if j == 1
                ylabel(['cell ' num2str(this_cell)])
            end
        end
%         g_vals(i) = l23_cell_analyzed17(this_cell).g;
    else
        for j = 1:1%length(l23_cell_analyzed16(this_cell).spike_data)
        

                plot([10 25 50 100 150],l23_cell_analyzed16(this_cell).spike_data(j).num_spike_means)

            hold on
        end
    end
    
end


%% spike timing
figure;
% cell_select = 26:38;
% cell_select = setdiff(cell_select,[29 31 32]);
% colors = parula(length(cell_select));
count = 1;
for i = 1:length(cell_select)
    this_cell = cell_select(i);
%     subplot(ceil(length(cell_select)/4),4,i)
    if l23_cells(this_cell).do_cc
        for j = 1:length(l23_cell_analyzed_fulldata_longer(this_cell).voltage_data)
            subplot(length(cell_select),10,(i-1)*10 + j)
            plot(l23_cell_analyzed_fulldata_longer(this_cell).voltage_data(j).powers(1:size(l23_cell_analyzed_fulldata_longer(this_cell).glm_sim.spike_time_means,2)),...
                l23_cell_analyzed_fulldata_longer(this_cell).voltage_data(j).spike_times_means(1:size(l23_cell_analyzed_fulldata_longer(this_cell).glm_sim.spike_time_means,2))/20,'k.-')

            hold on
            plot(l23_cell_analyzed_fulldata_longer(this_cell).voltage_data(j).powers(1:size(l23_cell_analyzed_fulldata_longer(this_cell).glm_sim.spike_time_means,2)),...
                l23_cell_analyzed_fulldata_longer(this_cell).glm_sim.spike_time_means(j,:)/20,'bo-')
%             hold on
%              plot(l23_cell_analyzed21(this_cell).voltage_data(j).powers(1:size(l23_cell_analyzed21(this_cell).glm_sim.spike_count_means,2)),...
%                 l23_cell_analyzed21(this_cell).glm_sim.spike_count_means(j,:),'m')
            
            count = count + 1;
            title(mat2str(l23_cell_analyzed_fulldata_longer(this_cell).voltage_data(j).location))
            ylim([0 400]/20)
            xlim([0 110])
            if j == 1
                ylabel(['cell ' num2str(this_cell)])
            end
        end
%         g_vals(i) = l23_cell_analyzed17(this_cell).g;
    else
        for j = 1:1%length(l23_cell_analyzed16(this_cell).spike_data)
        

                plot([10 25 50 100 150],l23_cell_analyzed24(this_cell).spike_data(j).num_spike_means)

            hold on
        end
    end
    
end
%%
figure;

subplot(151)
plot([l23_cell_analyzed17(cell_select).g; l23_cell_analyzed18(cell_select).g])
title('g')
subplot(152)
plot([l23_cell_analyzed17(cell_select).v_th; l23_cell_analyzed18(cell_select).v_th])
title('v-th')
subplot(153)
plot([l23_cell_analyzed17(cell_select).v_reset; l23_cell_analyzed18(cell_select).v_reset])
title('v-reset')
subplot(154)
plot([l23_cell_analyzed17(cell_select).gain; l23_cell_analyzed18(cell_select).gain])
title('gain')
subplot(155)
plot([l23_cell_analyzed17(cell_select).th_gain_ratio; l23_cell_analyzed18(cell_select).th_gain_ratio])
title('rheobase')
ylim([0 5000])
        
%% predict spikes from spatial model

% compute mean shape
analyzed_cells = 1;
cells_for_shape = find([l23_cells.do_vc]);
z_slices = size(analyzed_cells(1).current_data.shape_max,3);
average_shape = zeros(9,9,7);
count = 1;
for i = 1:length(cells_for_shape)
    cell_i = cells_for_shape(i);
    this_shape = analyzed_cells(cell_i).current_data.shape_max;
    this_shape = this_shape/max(this_shape(:));
    average_shape = average_shape + this_shape;
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
% get shape or mean shape

%% get voltage shape from cell 26

voltage_data = l23_cell_analyzed17(26).voltage_data;

no_spike_trials = [];
for i = 1:length(voltage_data)
    for j = 1:length(voltage_data(i).num_spike_means)
        if voltage_data(i).num_spike_means(j) == 0
            if isempty(no_spike_trials)
                no_spike_trials = voltage_data(i).data{j};
            else
                no_spike_trials = [no_spike_trials; voltage_data(i).data{j}];
            end
        end
    end
end

[u,s,v] = svd(bsxfun(@minus,no_spike_trials,no_spike_trials(:,1)));
if mean(v(:,1)) < 0
    current_shape = -v(:,1);
else
    current_shape = v(:,1);
end

current_shape = current_shape - current_shape(1);
ccurrent_shape = current_shape/max(current_shape);
current_shape_sv = s(1)/sum(diag(s));




figure; plot(current_shape)

%% check spikes

% cell_select = 26:38;
cell_select = 33:35
for i = 1:length(cell_select)
    
    ii = cell_select(i);
    v_data = l23_cell_analyzed20(ii).voltage_data;
    for j = 1:length(v_data)
        h = figure;
        plot_trace_stack_grid(v_data(j).data,Inf,1,0,[],[],[],v_data(j).spike_times);
        title(['cell ' num2str(ii) ', loc: ' mat2str(v_data(j).location)])
        set(gcf,'Position',[675 0 500 975])
    end
    waitfor(h)
end








