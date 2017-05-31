clear l23_cells
base_cell.start_trial = 1;
base_cell.cc_spike_thresh = 8;
base_cell.ca_num_spike_locs = 5;
base_cell.ca_spike_thresh = 6;
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
base_cell.fit_cc = 0;
base_cell.glm_sim_scale = 1;
base_cell.trial_dur = .010*20000;
base_cell.use_shape = 1;
base_cell.fit_locs = [];

this_cell = 1;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice1_cell1.mat';
l23_cells(this_cell).exclude_trials = [1 2 3];
l23_cells(this_cell).fluor = 4488-3768;


this_cell = 2; % not great location choices
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice1_cell2.mat';
% l23_cells(this_cell).cc_spike_thresh = 10;
l23_cells(this_cell).exclude_trials = [1 4];
l23_cells(this_cell).fluor = 2548-2284;


this_cell = 3;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice2_cell1.mat';
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).ca_spike_thresh = 5.5;
l23_cells(this_cell).exclude_trials = [1 3];
l23_cells(this_cell).fluor = 3256-2504;

this_cell = 4;
l23_cells(this_cell) = base_cell;% really only looking at first spikes here
l23_cells(this_cell).filename = '4_2_slice2_cell2.mat';
l23_cells(this_cell).ca_spike_thresh = 10;
l23_cells(this_cell).exclude_trials = [1 2 3 4 5];
l23_cells(this_cell).fluor = 2712-2332;

this_cell = 5; % really only looking at first spikes here
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice2_cell3.mat';
l23_cells(this_cell).ca_spike_thresh = 8;
l23_cells(this_cell).exclude_trials = [];
l23_cells(this_cell).fluor = 2797-2366;

this_cell = 6; % doublets
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice3_cell1.mat';
l23_cells(this_cell).exclude_trials = [1 4 5];
l23_cells(this_cell).fluor = 2653-2437;
l23_cells(this_cell).ca_spike_thresh = 10;

this_cell = 7; % doublets - hard to detect second spikes on this cell
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice3_cell2.mat';
l23_cells(this_cell).ca_spike_thresh = 15;
% base_cell.exclude_trials = [1 4 5];
l23_cells(this_cell).fluor = 2763-2331;

this_cell = 8; % very quick second spikes!
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_2_slice4_cell1.mat';
l23_cells(this_cell).ca_spike_thresh = 8;
l23_cells(this_cell).fluor = 2804-2411;

this_cell = 9; % NO SPIKES - CONTROL CELL
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice4_cell2.mat';
l23_cells(this_cell).ca_spike_thresh = 50; % no spikes
l23_cells(this_cell).do_vc = 0;
% l23_cells(this_cell).exclude_trials = [1 2 3 4 5];
l23_cells(this_cell).fluor = 2544-2528;

this_cell = 10;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice1_cell1.mat';
l23_cells(this_cell).start_trial = 3;


this_cell = 11; %doublets - maybe good for testing doublets
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice1_cell2.mat';
l23_cells(this_cell).start_trial = 3;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).ca_spike_thresh = 8;

this_cell = 12; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice1_cell3.mat';
l23_cells(this_cell).start_trial = 2;
l23_cells(this_cell).do_vc = 0;

this_cell = 13; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice2_cell1.mat';
l23_cells(this_cell).start_trial = 3;
l23_cells(this_cell).ca_spike_thresh = 10;

this_cell = 14; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice2_cell2.mat';
l23_cells(this_cell).ca_spike_thresh = 10;

this_cell = 15; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice3_cell1.mat';
l23_cells(this_cell).ca_spike_thresh = 20;

this_cell = 16;  % not a good cell for lif-glm
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice4_cell1.mat';
l23_cells(this_cell).ca_spike_thresh = 7;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).ca_spike_thresh = 15;

this_cell = 17; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice4_cell2next.mat';
l23_cells(this_cell).ca_spike_thresh = 10;

this_cell = 18; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_3_slice5_cell1next.mat';

this_cell = 19; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice1_cell1next.mat';
l23_cells(this_cell).ca_spike_thresh = 20;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).start_trial = 3;
l23_cells(this_cell).fluor = 3004-2788;

this_cell = 20; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice1_cell2.mat';
% l23_cells(this_cell).ca_spike_thresh = 15;
l23_cells(this_cell).start_trial = 2;
l23_cells(this_cell).fluor = 292-189;

this_cell = 21; % NOT A LOT OF SPIKES
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice1_cell3.mat';
% l23_cells(this_cell).ca_spike_thresh = 8;
l23_cells(this_cell).start_trial = 2;
l23_cells(this_cell).fluor = 261-202;

this_cell = 22; % START FLUOR HERE??
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice2_cell1next.mat';
% l23_cells(this_cell).ca_spike_thresh = 8;
l23_cells(this_cell).start_trial = 3;
l23_cells(this_cell).fluor = 590-220;

this_cell = 23; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice2_cell2.mat';
% l23_cells(this_cell).ca_spike_thresh = 15;
l23_cells(this_cell).start_trial = 2;
l23_cells(this_cell).fluor = 720-223;

this_cell = 24; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice3_cell1next.mat';
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).fluor = 878-197;
l23_cells(this_cell).ca_spike_thresh = 12;

this_cell = 25; 
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_6_slice3_cell2next.mat';
l23_cells(this_cell).start_trial = 2;
l23_cells(this_cell).fluor = 857-217;
l23_cells(this_cell).ca_spike_thresh = 10;

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
l23_cells(this_cell).cc_spike_thresh = 20;

this_cell = 28;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '3_30_slice2_cell2.mat';
l23_cells(this_cell).ca_num_spike_locs = 10;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 10;
l23_cells(this_cell).cc_spike_thresh = 20;

this_cell = 29;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '3_30_slice2_cell3.mat';
l23_cells(this_cell).ca_num_spike_locs = 10;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 10;
l23_cells(this_cell).cc_spike_thresh = 20;

this_cell = 30;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '3_30_slice4_cell2.mat';
l23_cells(this_cell).ca_num_spike_locs = 10;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 10;
l23_cells(this_cell).cc_spike_thresh = 20;

this_cell = 31;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_21_slice3_cell2.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 10;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).cc_spike_thresh = 20;
l23_cells(this_cell).ca_spike_thresh = -15;

this_cell = 32;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_25_slice1_cell1.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 6;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).cc_spike_thresh = 20;
l23_cells(this_cell).ca_spike_thresh = -25;

this_cell = 33;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_25_slice1_cell2.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 6;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).cc_spike_thresh = 20;
l23_cells(this_cell).ca_spike_thresh = -25;

this_cell = 34;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_25_slice2_cell1.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 6;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).cc_spike_thresh = 20;

this_cell = 35;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_25_slice2_cell2.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 6;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).cc_spike_thresh = 20;

this_cell = 36;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_25_slice2_cell3.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 6;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).cc_spike_thresh = 20;

this_cell = 37;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_25_slice3_cell1.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 6;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).cc_spike_thresh = 20;

this_cell = 38;
l23_cells(this_cell) = base_cell;
l23_cells(this_cell).filename = '4_25_slice3_cell2.mat';
l23_cells(this_cell).ca_num_spike_locs = 1;
l23_cells(this_cell).do_cc = 1;
l23_cells(this_cell).cc_num_spike_locs = 6;
l23_cells(this_cell).do_vc = 0;
l23_cells(this_cell).cc_spike_thresh = 20;



%%

cell_to_run = 1:38;%cell_select;%[26:38];
cell_to_run = setdiff(cell_to_run,[16]);
% cell_to_run = 26:30;
cell_to_run = 27;
cell_to_run = find([l23_cells.do_vc]);
% cell_analyzed_bu = cell_analyzed;
% clear l23_cell_analyzed13 <-- multispike, 14 is one spike
% clear l23_cell_analyzed_10ms_fulldata_noshape
for i = 1:length(cell_to_run)
    
    this_cell = cell_to_run(i)
    l23_cell_analyzed_10ms_fulldata_centavgshape(this_cell) = ...
        analyze_opto_response(l23_cells(this_cell));


    
end


%% l5 cells
base_cell.start_trial = 1;
base_cell.cc_spike_thresh = 8;
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
l5_cells(this_cell).cc_spike_thresh = 16;

this_cell = 2; 
l5_cells(this_cell) = base_cell;
l5_cells(this_cell).filename = '4_3_slice4_cell3.mat';

this_cell = 3; 
l5_cells(this_cell) = base_cell;
l5_cells(this_cell).filename = '4_3_slice5_cell2.mat';

this_cell = 4; 
l5_cells(this_cell) = base_cell;
l5_cells(this_cell).filename = '4_6_slice2_cell3.mat';
l5_cells(this_cell).cc_spike_thresh = 12;

this_cell = 5; 
l5_cells(this_cell) = base_cell;
l5_cells(this_cell).filename = '4_6_slice3_cell3next.mat';
l5_cells(this_cell).hpass_filt = 1;
% l5_cells(this_cell).cc_spike_thresh = 12;
l5_cells(this_cell).do_vc = 0;

this_cell = 6; 
l5_cells(this_cell) = base_cell;
l5_cells(this_cell).filename = '4_6_slice3_cell4next.mat';
l5_cells(this_cell).hpass_filt = 1;
l5_cells(this_cell).cc_spike_thresh = 12;



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

cells_to_plot = setdiff(cells_to_plot,[5 27:30]);
% cells_to_plot = 25;
this_analysis = l23_cell_analyzed_10ms_fulldata_noshape_constg;
% cells_to_plot = 27:30;
z_slices = size(this_analysis(cells_to_plot(1)).current_data.shape_max,3);
% l23_average_shape = zeros(9,9,z_slices);

z_depths = {'-90','-50','-20','0','20','50','90'};
% z_depths = {'-60','-40','-20', '-10', '0', '10', '20','40','60'};
count = 1;
for i = 1:length(cells_to_plot)
    cell_i = cells_to_plot(i);
    this_shape = this_analysis(cell_i).current_data.shape_max;
    this_shape = this_shape/max(this_shape(:));
%     l23_average_shape = l23_average_shape + this_shape;

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

% l23_average_shape = l23_average_shape/max(l23_average_shape(:));

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
cells_to_plot = setdiff(cells_to_plot,[5 6 26 27 29 30 ]);
avg_template = [];
colors = parula(100);
svals = zeros(length(cells_to_plot),1);
for i = 1:length(cells_to_plot)
    ii = cells_to_plot(i);
    svals(i) = l23_cell_analyzed_10ms_fulldata_noshape_constg(ii).current_data.current_shape_sv;
end
% svals = floor(svals/max(svals)*100);
svals = floor(svals*100);
for i = 1:length(cells_to_plot)
    ii = cells_to_plot(i);
    if i == 1
        avg_template = l23_cell_analyzed_10ms_fulldata_noshape_constg(ii).current_data.current_shape;
    else
        avg_template = avg_template + l23_cell_analyzed_10ms_fulldata_noshape_constg(ii).current_data.current_shape;
    end
%     subplot(length(cells_to_plot),1,i)
    plot((1:length(l23_cell_analyzed_10ms_fulldata_noshape_constg(ii).current_data.current_shape))/20000,...
        -l23_cell_analyzed_10ms_fulldata_noshape_constg(ii).current_data.current_shape, 'Color',[.33 .33 .33])%colors(svals(i),:))
    hold on
end
avg_template = avg_template/length(cells_to_plot);
plot((1:length(avg_template))/20000,-avg_template,'k','linewidth',3)
hold off
xlabel('time (sec)')
ylabel('current (a.u.)')
title('est. chrome current shape - 3 msec pulse - l23pyr (N = 17)')
xlim([0 .04])
axis off


%%

% analyzed_cells = these_cells_analyzed;
cell_select = 26:38;

opt_rheobase = zeros(size(cell_select));
mean_time = zeros(size(cell_select));

for i = 1:length(cell_select)
    
    this_cell = cell_select(i);

    op_rheobase(i) = these_cells_analyzed(this_cell).th_gain_ratio;


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

%     semilogx(g,[l23_cell_analyzed5(this_cell).glm_out.dev these_cells_analyzed(this_cell).glm_out.dev l23_cell_analyzed6(this_cell).glm_out.dev],'color',colors(i,:))
plot(these_cells_analyzed(this_cell).glm_params.g,[these_cells_analyzed(this_cell).glm_out.dev],'color',colors(i,:))


%     semilogx(g,[l23_cell_analyzed5(this_cell).glm_out.dev l23_cell_analyzed_fulldata_longer(this_cell).glm_out.dev l23_cell_analyzed6(this_cell).glm_out.dev],'color',colors(i,:))
    hold on
    
end


%% num_spikes

figure;
% cell_select = 1:38;
cell_select = find(~[l23_cells.do_cc]);
% cell_select = setdiff(cell_select,[29 31 32]);
% cell_select = cell_select([1 2 5 6]);
colors = parula(length(cell_select));
count = 1;
these_cells_analyzed = l23_cell_analyzed_10ms_fulldata_noshape;
for i = 1:length(cell_select)
    this_cell = cell_select(i);
%     subplot(ceil(length(cell_select)/4),4,i)
    if l23_cells(this_cell).do_cc% && l23_cells(this_cell).do_vc
        spike_data = these_cells_analyzed(this_cell).voltage_data;
    else
        spike_data = these_cells_analyzed(this_cell).spike_data;
    end
    for j = 1:length(spike_data)
        subplot(length(cell_select),10,(i-1)*10 + j)
        plot(spike_data(j).powers(1:size(these_cells_analyzed(this_cell).glm_sim.spike_count_means,2)),...
            spike_data(j).num_spike_means(1:size(these_cells_analyzed(this_cell).glm_sim.spike_count_means,2)),'ko')

        hold on
        plot(spike_data(j).powers(1:size(these_cells_analyzed(this_cell).glm_sim.spike_count_means,2)),...
            these_cells_analyzed(this_cell).glm_sim.spike_count_means(j,:),'b.-')
%             hold on
%              plot(l23_cell_analyzed21(this_cell).voltage_data(j).powers(1:size(l23_cell_analyzed21(this_cell).glm_sim.spike_count_means,2)),...
%                 l23_cell_analyzed21(this_cell).glm_sim.spike_count_means(j,:),'m')

        count = count + 1;
%         title(mat2str(spike_data(j).location))
        ylim([0 3.0])
        xlim([0 110])
        set(gca,'xticklabel',{})
        set(gca,'yticklabel',{})
        if j == 1
            ylabel([num2str(this_cell)])
        end
    end
%         g_vals(i) = l23_cell_analyzed17(this_cell).g;
    
    
end


% spike timing
figure;
% cell_select = 26:38;
% cell_select = setdiff(cell_select,[29 31 32]);
% colors = parula(length(cell_select));
count = 1;
for i = 1:length(cell_select)
    this_cell = cell_select(i);
%     subplot(ceil(length(cell_select)/4),4,i)
    if l23_cells(this_cell).do_cc% && l23_cells(this_cell).do_vc
        spike_data = these_cells_analyzed(this_cell).voltage_data;
    else
        spike_data = these_cells_analyzed(this_cell).spike_data;
    end
    for j = 1:length(spike_data)
        subplot(length(cell_select),10,(i-1)*10 + j)
        plot(spike_data(j).powers(1:size(these_cells_analyzed(this_cell).glm_sim.spike_time_means,2)),...
            spike_data(j).spike_times_means(1:size(these_cells_analyzed(this_cell).glm_sim.spike_time_means,2))/20,'ko')
        hold on
        errorbar(spike_data(j).powers(1:size(these_cells_analyzed(this_cell).glm_sim.spike_time_means,2)),...
            these_cells_analyzed(this_cell).glm_sim.spike_time_means(j,:)/20,...
            these_cells_analyzed(this_cell).glm_sim.spike_time_std(j,:)/20,'b.-')
%             hold on
%              plot(l23_cell_analyzed21(this_cell).voltage_data(j).powers(1:size(l23_cell_analyzed21(this_cell).glm_sim.spike_count_means,2)),...
%                 l23_cell_analyzed21(this_cell).glm_sim.spike_count_means(j,:),'m')

        count = count + 1;
%         title(mat2str(spike_data(j).location))
        ylim([0 400]/20)
        xlim([0 110])
        set(gca,'xticklabel',{})
        set(gca,'yticklabel',{})
        if j == 1
            ylabel([num2str(this_cell)])
        end
    end
%         g_vals(i) = l23_cell_analyzed17(this_cell).g;

    
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
cell_select = 1:5;
for i = 1:length(cell_select)
    
    ii = cell_select(i);
    v_data = l23_cell_analyzed_preprocessonly_new(ii).spike_data;
    for j = 1:length(v_data)
        h = figure;
        plot_trace_stack_grid(v_data(j).data,Inf,1,0,[],[],[],v_data(j).spike_times);
        title(['cell ' num2str(ii) ', loc: ' mat2str(v_data(j).location)])
        xlim([0 .015])
        set(gcf,'Position',[675 0 500 975])
    end
    waitfor(h)
end


%% better shape template
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
%     this_shape = smoothn(this_shape,.1);
    [~, center_ind] = max(this_shape(:));
    center = zeros(1,3);
    [center(1), center(2), center(3)] = ind2sub([9,9,7],center_ind);
    offset = center - [5 5 4];
    this_shape_tmp = zeros(size(this_shape));
    for j = (1:9) - offset(1)
        if j > 0 && j < 10
        for k = (1:9) - offset(2)
            if k > 0 && k < 10
            for m = (1:7) - offset(3)
                if m > 0 && m < 8
                    this_shape_tmp(j,k,m) = this_shape(j+offset(1),...
                                                       k+offset(2),...
                                                       m+offset(3));
                end
            end
            end
        end
        end
    end
    this_shape = this_shape_tmp/max(this_shape_tmp(:));
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

gains_avg = zeros(9,9,7);
gains_count = zeros(9,9,7);
count = 1;
z_depths = [-20 0 20];
for i = 1:length(cell_select)
    this_cell = cell_select(i);
    
        this_gain = these_cells_analyzed(this_cell).gain;
        this_gain = this_gain/this_gain(1);
        for j = 1:length(spike_data)
            if this_gain(j) < 0
                this_gain(j) = 0;
            end
            this_loc = spike_data(j).location;
            this_loc(1:2) = this_loc(1:2)/10 + 5;
            this_loc(3) = find(round(this_loc(3),-1) == z_depths) + 2;
            gains_avg(this_loc(1),this_loc(2),this_loc(3)) = ...
                gains_count(this_loc(1),this_loc(2),this_loc(3)) + ...
                this_gain(j);
            gains_count(this_loc(1),this_loc(2),this_loc(3)) = ...
                gains_count(this_loc(1),this_loc(2),this_loc(3)) + 1;
        end
%         g_vals(i) = l23_cell_analyzed17(this_cell).g;
    
end
gains_avg = gains_avg./gains_count;
gains_avg_norm = gains_avg/max(gains_avg(:));
gains_avg_norm(gains_count(:) < 2) = NaN;

%% plot shapes
figure
z_slices = size(shapes,3);
for i = 1:size(shapes,4);
for j = 1:z_slices
    subplot(z_slices,size(shapes,4),(j-1)*(size(shapes,4)) + i)
    imagesc(shapes(:,:,j,i))
%     axis off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    caxis([0 1])
    
    if j == 1
        title(shape_names{i})
    end
    ylabel(['z = ' z_depths{j}])
end
end


%% interpolated shapes

figure
cells_to_plot = find([l23_cells.do_vc]);
cells_to_plot = setdiff(cells_to_plot,[5]);
this_analysis = l23_cell_analyzed_10ms_preprocess_only_raw_shape;
% cells_to_plot = 27:30;
% z_slices = size(this_analysis(cells_to_plot(1)).current_data.shape_max,3);
l23_average_shape = zeros(81,81,181);
z_depths = [-40 -30 -20 -10 0 10 20 30 40];
z_slices = length(z_depths);
all_shapes = nan(numel(l23_average_shape),length(cells_to_plot));
% z_depths = {'-60','-40','-20', '-10', '0', '10', '20','40','60'};
count = 1;
for i = 1:length(cells_to_plot)
    cell_i = cells_to_plot(i);
    this_shape = this_analysis(cell_i).current_data.upres_shape;
%     this_shape = this_shape/max(this_shape(:));
    this_shape = padarray(this_shape,[0 0 (181-size(this_shape,3))/2],NaN);
%     l23_average_shape = l23_average_shape + this_shape;
    
    % centered shape
     [~, center_ind] = max(this_shape(:));
    center = zeros(1,3);
    [center(1), center(2), center(3)] = ind2sub([81,81,181],center_ind);
    offset = center - [41 41 91];
    this_shape_tmp = nan(size(this_shape));
    for j = (1:81) - offset(1)
        if j > 0 && j < 82
        for k = (1:81) - offset(2)
            if k > 0 && k < 82
            for m = (1:181) - offset(3)
                if m > 0 && m < 182
                    this_shape_tmp(j,k,m) = nansum([this_shape_tmp(j,k,m) ...
                                                       this_shape(j+offset(1),...
                                                       k+offset(2),...
                                                       m+offset(3))]);
                end
            end
            end
        end
        end
    end
    this_shape = this_shape_tmp;%/max(this_shape_tmp(:));
    all_shapes(:,i) = this_shape(:);
    
    z_inds = ceil((z_depths - ...
        -90)/current_data.upres);
    l23_plot_shape = this_shape(:,:,z_inds);
    for j = 1:z_slices
        subplot(z_slices,length(cells_to_plot)+1,(j-1)*(length(cells_to_plot)+1) + i + 1)
        imagesc(l23_plot_shape(:,:,j))
        caxis([0 1200])
        axis off
%         axis image
        count = count + 1;
        if j == 1
            title(['cell ' num2str(cell_i)])
        end
    end
end
l23_average_shape = ...
        reshape(nanmean(all_shapes,2),[81 81 181]);
    

% l23_average_shape = l23_average_shape/length(cells_to_plot);%max(l23_average_shape(:));
l23_average_shape_norm = l23_average_shape/max(l23_average_shape(:));
figure
for j = 1:z_slices
    subplot(z_slices,length(cells_to_plot)+1,(j-1)*(length(cells_to_plot)+1) + 1)
    imagesc(l23_average_shape(:,:,z_inds(j)))
%     axis off
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    caxis([0 1])
    axis off
%     axis image
    if j == 1
        title('mean cell')
    end
    ylabel(['z = ' num2str(z_depths(j))])
end

% l23_average_shape_raw = l23_average_shape/max(l23_average_shape(:));



%% compare gains vs. shape

cells_to_plot = find([l23_cells.do_vc]);
cells_to_plot = setdiff(cells_to_plot,[21 26:38]);
gains = [];
shape_scale = [];
cell_id = [];

fulldata_struct = l23_cell_analyzed_10ms_fulldata_centavgshape_constg;
upres_struct = l23_cell_analyzed_10ms_fulldata_centavgshape_constg;
figure
colors = hsv(ceil(length(cells_to_plot)*1.75));
avg_shape = l23_average_shape_norm;
upres_x_locs = -40:1:40;
upres_y_locs = -40:1:40;
upres_z_locs = -90:1:90;
for i = 1:length(cells_to_plot)
%     h = figure
    cell_i = cells_to_plot(i);
    this_cell = fulldata_struct(cell_i);
    this_gain = this_cell.gain;
    
    if this_cell.do_cc
        these_spikes = this_cell.voltage_data;
    else
        these_spikes = this_cell.spike_data;
    end
    
    this_shape = upres_struct(cell_i).current_data.upres_shape;
    this_shape = this_shape/max(this_shape(:));
    
    
%     this_shape = l23_average_shape_norm;
    location_gains = zeros(size(this_gain));
%     upres_x_locs = upres_struct(cell_i).current_data.upres_x_locs;
%     upres_y_locs = upres_struct(cell_i).current_data.upres_y_locs;
%     upres_z_locs = upres_struct(cell_i).current_data.upres_z_locs;
    for k = 1:length(these_spikes)

        this_loc = these_spikes(k).location;

        location_ind = [find(this_loc(1) == upres_x_locs) ...
                        find(this_loc(2) == upres_y_locs) ...
                        find(this_loc(3) == upres_z_locs)];
                    
        location_gains(k) = ...
            this_shape(location_ind(1),location_ind(2),location_ind(3))*this_gain(1);
        avg_loc_gain(k) = ...
            avg_shape(location_ind(1),location_ind(2),location_ind(3));
    end
    
    gains = [gains; this_gain(1)*avg_loc_gain'];
    shape_scale = [shape_scale; location_gains];
    cell_id = [cell_id; cell_i*ones(size(this_gain))];
    
    lg_bu = location_gains;
    gain_bu = this_gain;
    gain_bu(gain_bu < 0) = 0;
    location_gains(this_gain < 0) = [];
    this_gain(this_gain < 0) = [];
%     this_gain(location_gains >= 1) = [];
%     location_gains(location_gains >= 1) = [];
    if length(location_gains) > 2
%         h = figure;
        scatter(location_gains,this_gain,[],repmat(colors(i,:),length(location_gains),1))
        hold on
%         lsline
        p = polyfit(location_gains,this_gain,1);
        plot(location_gains,p(1)*location_gains + p(2),'color',colors(i,:))
        hold on
        scatter(location_gains,this_gain,[],repmat(colors(i,:),length(location_gains),1))
        xlim([0 1.2]*1)
        ylim([0 .02])
        
        scatter(lg_bu,gain_bu)
        hold on
%         waitfor(h)
    end
%     waitfor(h)
end

hold off
gains_regres = gains;
shape_scale_regres = shape_scale;

gains(gains < 0) = 0;
% shape_scale(shape_scale < 0) = 0;

shape_scale_regres(gains_regres < 0) = [];
gains_regres(gains_regres < 0) = [];
% gains_regres(shape_scale_regres >= 1) = [];
% shape_scale_regres(shape_scale_regres >= 1) = [];

figure
gscatter(shape_scale,gains,cell_id)
hold on
plot(0:.001:0.02,(0:.001:0.02))
% xlim([0 1.1])
% ylim([0 1.1])

ylabel('lif-glm location gain (avg. shape)')
xlabel('v-clamp shape gain')

%% ks test on neurons

analysis_to_run = {l23_cell_analyzed_10ms_fulldata_noshape_fitvth,...
    l23_cell_analyzed_10ms_fulldata_noshape,...
    l23_cell_analyzed_10ms_fulldata_noshape_fitvth_constg,...
    l23_cell_analyzed_10ms_fulldata_noshape_constg,...
    l23_cell_analyzed_10ms_fulldata_centavgshape_constg_fitvth,...
    l23_cell_analyzed_10ms_fulldata_centavgshape_constg};
% this_analysis = l23_cell_analyzed_10ms_fulldata_centavgshape_constg_fitvth;
glm_type = {'full model','const g','const vth','const g and vth','template shape const g',...
    'template shape const g const vth'};
% cells_to_do = setdiff(1:38,[16]);
cells_to_do = find(~[l23_cells.do_cc]);
figure
the_handles = [];
colors = lines(length(analysis_to_run));
for jj = 1: length(analysis_to_run)
    this_analysis = analysis_to_run{jj};
    all_z = [];
    for ii = 1:length(cells_to_do);

        cell_i = cells_to_do(ii);
        this_lambda = this_analysis(cell_i).lambda';
        responses = zeros(size(this_lambda));
    %     this_lambda = this_lambda(:);
        if this_analysis(cell_i).do_cc
            num_spike_locs = this_analysis(cell_i).cc_num_spike_locs;

            these_spikes = this_analysis(cell_i).voltage_data;
        else
            num_spike_locs = this_analysis(cell_i).ca_num_spike_locs;
            these_spikes = this_analysis(cell_i).spike_data;
        end
        count = 1;
        for k = 1:num_spike_locs

    %         k = fit_locs(kk);
            data_spike_times = these_spikes(k).spike_times;
            powers = these_spikes(k).powers;
            powers = powers(1:end-1);%-1 % don't do highest power - never useful

            for i = 1:length(powers)


                for j = 1:length(data_spike_times{i})
    %                 responses(:,count) = zeros(1,length(current_template(1:1:end)));

                    if ~isempty(data_spike_times{i}{j})
                        responses(floor(data_spike_times{i}{j}/1),count) = 1;

                    end
                    count = count + 1;
                end
            end    
        end

    %     responses = responses(:);

        this_lambda_cum = cumsum(this_lambda,1);
        tau = [];
        for j = 1:size(responses,2);
            this_tau = diff(this_lambda_cum([1 find(responses(:,j)')],j))';
            tau = [tau this_tau];
        end
        z = 1 - exp(-tau);
        all_z = [all_z z];
    %     if ~isempty(z)
    %         return
    %     end

    end

    [z_trans_sorted z_sort_ind] = sort(all_z);

    b = ((1:length(all_z)) - 0.5)/length(all_z);

%     figure
    h = plot(z_trans_sorted,b,'color',colors(jj,:));
    the_handles = [the_handles h];
    hold on;
%     plot(z_trans_sorted,b-1.36/sqrt(length(all_z)),'--',...
%                                     z_trans_sorted,b+1.36/sqrt(length(all_z)),'--','Color',colors(jj,:))
end
hold on
x = 0:.1:1;
plot(x,x,'r')
title(['KS Test For LIF-GLM Version w/ bounds'])
% legend('Emperical CDF','CI Bound','CI Bound','Uniform CDF')
legend(the_handles,glm_type); 



%% plot intensity vs spiking stats for a single cell

cell_ind = 27;
cell_select = find(~[l23_cells.do_cc]);
load('l23_centered_upres_cell.mat')
upres_x_locs = -40:1:40;
upres_y_locs = -40:1:40;
upres_z_locs = -90:1:90;
this_cell_shape = l23_average_shape_norm;
figure
x_counts = [];
all_data_y_counts = [];
all_preds_y_counts = [];

x_times = [];
all_data_y_times = [];
all_preds_y_times = [];

all_data_y_jitter = [];
all_preds_y_jitter = [];

do_sim = 0;
for i = 1:length(cell_select)
    cell_ind = cell_select(i);

    this_cell_spike_data = l23_cell_analyzed_10ms_fulldata_centavgshape_constg(cell_ind).spike_data;
    this_cell_sim_data = l23_cell_analyzed_10ms_fulldata_centavgshape_constg(cell_ind).glm_sim;
%     this_cell_shape = l23_cell_analyzed_10ms_fulldata_centavgshape_constg(cell_ind).current_data.shape_svd;
%     this_cell_shape = this_cell_shape/max(this_cell_shape(:));
    
%     upres_x_locs = -40:10:40;
%     upres_y_locs = -40:10:40;
%     upres_z_locs = [-90 -50 -20 0 20 50 90];
    powers = [10 25 50 100];
    
    


    for k = 1:length(this_cell_spike_data)

        this_loc = round(this_cell_spike_data(k).location,-1);

        location_ind = [find(this_loc(1) == upres_x_locs) ...
                        find(this_loc(2) == upres_y_locs) ...
                        find(this_loc(3) == upres_z_locs)];

        if length(location_ind) < 3
            continue
        end
        shape_gain = ...
            this_cell_shape(location_ind(1),location_ind(2),location_ind(3));%*...
            %l23_cell_analyzed_10ms_fulldata_centavgshape_constg(cell_ind).gain(1);
        stim_intensity = shape_gain*powers;
        stim_intensity(stim_intensity < 0) = 0;
    %     for i = 1:length(powers)
        subplot(311)
        scatter(stim_intensity,this_cell_spike_data(k).num_spike_means(1:4),'b','filled', 'jitter','on', 'jitterAmount',0.5);
        good_inds = this_cell_spike_data(k).num_spike_means(1:4) > .4;
        hold on
        subplot(312)
        scatter(stim_intensity(good_inds),this_cell_spike_data(k).spike_times_means(good_inds)/20,'b','filled', 'jitter','on', 'jitterAmount',0.5);
        hold on
        subplot(313)
        scatter(stim_intensity(good_inds),this_cell_spike_data(k).spike_times_std(good_inds)/20,'b','filled', 'jitter','on', 'jitterAmount',0.5);
        hold on
        
        x_counts = [x_counts stim_intensity];
        all_data_y_counts = [all_data_y_counts this_cell_spike_data(k).num_spike_means(1:4)];
        x_times = [x_times stim_intensity(good_inds)];
        all_data_y_times = [all_data_y_times this_cell_spike_data(k).spike_times_means(good_inds)/20];
        all_data_y_jitter = [all_data_y_jitter this_cell_spike_data(k).spike_times_std(good_inds)/20];

        if do_sim
            subplot(311)
            scatter(stim_intensity,this_cell_sim_data.spike_count_means(k,:),'r','filled', 'jitter','on', 'jitterAmount',0.5);
            hold on
    %         good_inds = this_cell_sim_data.spike_count_means(k,:) > .1;
            subplot(312)
            scatter(stim_intensity(good_inds),this_cell_sim_data.spike_time_means(k,good_inds)/20,'r','filled', 'jitter','on', 'jitterAmount',0.5);
            hold on
            subplot(313)
            scatter(stim_intensity(good_inds),this_cell_sim_data.spike_time_std(k,good_inds)/20,'r','filled', 'jitter','on', 'jitterAmount',0.5);
            hold on
            all_preds_y_counts = [all_preds_y_counts this_cell_sim_data.spike_count_means(k,:)];
            all_preds_y_times = [all_preds_y_times this_cell_sim_data.spike_time_means(k,good_inds)/20];
            all_preds_y_jitter = [all_preds_y_jitter this_cell_sim_data.spike_time_std(k,good_inds)/20];
        end

        
        


    end

end

x_max = max(union(x_counts,x_times));
% figure
subplot(311)
% f = fit(x_counts',all_data_y_counts','smoothingspline','SmoothingParam',0.07);
xval = min(x_counts):.01:max(x_counts);
fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)));
[param]=sigm_fit(x_counts,all_data_y_counts,[0 1 NaN NaN],[],0);
plot(xval,fsigm(param,xval),'b','Linewidth',2)
hold on
% f = fit(x_counts',all_preds_y_counts','smoothingspline','SmoothingParam',0.07);
if do_sim
    xval = min(x_counts):.01:max(x_counts);
    fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)));
    [param]=sigm_fit(x_counts,all_preds_y_counts,[0 1 NaN NaN],[],0);
    plot(xval,fsigm(param,xval),'r','Linewidth',2)
end
set(gca,'xticklabels',{})
ylim([0 1.1])
xlim([0 x_max])
ylabel('Mean Spike Count')
title('Spiking Statistics')
subplot(312)

f = fit(x_times(~isnan(all_data_y_times))',all_data_y_times(~isnan(all_data_y_times))','exp1');%,'SmoothingParam',0.07);
h = plot(f,'b')
set(h,'linewidth',2)
hold on
if do_sim
    f = fit(x_times(~isnan(all_preds_y_times))',all_preds_y_times(~isnan(all_preds_y_times))','exp1');%,'SmoothingParam',0.07);
    h = plot(f,'r')
    set(h,'linewidth',2)
end
set(gca,'xticklabels',{})
ylim([0 10])
xlim([0 x_max])
ylabel('Mean First Spike Time (msec)')
subplot(313)
f = fit(x_times(~isnan(all_data_y_jitter))',all_data_y_jitter(~isnan(all_data_y_jitter))','exp1');%,'smoothingspline','SmoothingParam',0.07);
h = plot(f,'b')
set(h,'linewidth',2)
hold on
if do_sim
    f = fit(x_times(~isnan(all_preds_y_jitter))',all_preds_y_jitter(~isnan(all_preds_y_jitter))','exp1');%,'smoothingspline','SmoothingParam',0.07);
    h = plot(f,'r')
    set(h,'linewidth',2)
end
xlabel('Est. Stim Intensity (a.u., avg. shape gain x laser power x cell gain)')
ylim([0 2])
xlim([0 x_max])
ylabel('First Spike Time Std. Dev. (msec)')













