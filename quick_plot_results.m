
close all

%%

iter_id = trials;



workloc = 1; % 0 is lab analysis computer, 1 is macbook pro
switch workloc
    case 0
        mappingroot = '/media/shababo/data/';
    case 1
        mappingroot = '~/projects/mapping/data/';
end


trials = 1:2;
duration = .030;
loc_trials = {1:2}; % each entry is a vector of trials that belong to that entry index location (e.g. {3:9,10:17,18:26})

loc_ids = [];
for i = 1:length(loc_trials)
    loc_ids = [loc_ids i*ones(size(loc_trials{i}))];
end
map_id = exp_data.params.map_id;

group_colors = [0 0 1;
                1 .8 .8;
                .8 1 .8;
                1 0 0;
                0 1 0];
            
group_colors(:) = 0;
            
group_names = {'undefined_cells','potentially_disconnected_cells',...
                'potentially_connected_cells','dead_cells','alive_cells'};


%%
all_locs = unique(loc_ids);
for i = length(all_locs)
    loc_id = all_locs(i);
    cell_group_list = exp_data.cells_targets.cell_group_list{loc_id};
    cell_group_locs = exp_data.cells_targets.cell_locations(cell_group_list,:);
    multispot_targs = exp_data.cells_targets.target_locations_selected{loc_id};
    single_targs = exp_data.cells_targets.target_locations_nuclei{loc_id};

    figure; subplot(121)
    scatter(cell_group_locs(:,2),-cell_group_locs(:,1),20*ones(size(cell_group_locs(:,1))),'filled'); hold on
    scatter(multispot_targs(:,2),-multispot_targs(:,1),8*ones(size(multispot_targs(:,1))),'r','filled')

    subplot(122)
    scatter(cell_group_locs(:,2),-cell_group_locs(:,1),20*ones(size(cell_group_locs(:,1))),'filled'); hold on
    scatter(single_targs(:,2),-single_targs(:,1),8*ones(size(single_targs(:,1))),'r','filled')
end

%%



iter = length(exp_data.design.(group_names{1}){loc_id});
figure
for i = 1:25
    subplot(3,3,i)
     iter = i
    for i = 1:length(group_names) 
    
        
       
    these_cells = logical(exp_data.design.(group_names{i}){loc_id}{iter});
    scatter(cell_group_locs(these_cells,2),-cell_group_locs(these_cells,1),...
        20*ones(size(cell_group_locs(these_cells,1))),...
        repmat(group_colors(i,:),sum(these_cells),1),...
        'filled');
    hold on
    end
end

%% plot gamma path


iter = length(exp_data.design.(group_names{1}){loc_id});
figure
for i = 1:25
    subplot(3,3,i)
     iter = i
    for i = 1:length(group_names) 
    
        
       
    these_cells = logical(exp_data.design.(group_names{i}){loc_id}{iter});
    scatter(cell_group_locs(these_cells,2),-cell_group_locs(these_cells,1),...
        real(exp_data.design.gamma_path{loc_id}(these_cells,iter))*20 + eps,...
        repmat(group_colors(i,:),sum(these_cells),1),...
        'filled');
    hold on
    end
end

%% plot gain path


iter = length(exp_data.design.(group_names{1}){loc_id});
figure
for i = 1:25
    subplot(3,3,i)
     iter = i
    for i = 1:length(group_names) 
    
        
       
    these_cells = logical(exp_data.design.(group_names{i}){loc_id}{iter});
    scatter(cell_group_locs(these_cells,2),-cell_group_locs(these_cells,1),...
        real(exp_data.design.gain_path{loc_id}(these_cells,iter))*1000 + eps,...
        repmat(group_colors(i,:),sum(these_cells),1),...
        'filled');
    hold on
    end
end


%% plot nuclei

% nuc_locs_img = nuc_locs(:,1:3);
nuc_locs_img = data.trial_metadata(end).nuclear_locs(:,1:3);
nuc_locs_img(:,1:2) = nuc_locs_img(:,1:2)/exp_data.image_um_per_px;
nuc_locs_img(:,3) = nuc_locs_img(:,3)/exp_data.stack_um_per_slice;
nuc_locs_img = bsxfun(@plus,nuc_locs_img,[exp_data.image_zero_order_coord' 0]);
nuc_locs_img(:,1:2) = nuc_locs_img(:,[2 1]);


plot_nuclear_detect_3D([mappingroot exp_data.params.map_id '_stack.tif'],nuc_locs_img');

%% plot targets

nuc_locs_img = [];

for i = 1:length(trials)
    
    this_trial = trials(i);
    these_targs = data.trial_metadata(this_trial).stim_key;
    for k = 1:size(these_targs,3)
        for j = 1:size(these_targs,1)
            if ~any(isnan(these_targs(j,:,k)))
                nuc_locs_img(end+1,:) = these_targs(j,:,k);
            end
        end
    end
    
    
end

nuc_locs_img(:,1:2) = nuc_locs_img(:,1:2)/exp_data.image_um_per_px;
nuc_locs_img(:,3) = nuc_locs_img(:,3)/exp_data.stack_um_per_slice;
nuc_locs_img = bsxfun(@plus,nuc_locs_img,[exp_data.image_zero_order_coord' 0]);
nuc_locs_img(:,1:2) = nuc_locs_img(:,[2 1]);


plot_nuclear_detect_3D([mappingroot exp_data.params.map_id '_stack.tif'],nuc_locs_img');


%%

this_seq = cell(length(trials),1);
this_stim_key = cell(length(trials),1);
power_curve_num = cell(length(trials),1);
stim_starts = cell(length(trials),1);

stims_per_trial = [];

full_stim_key = [];
clear full_seq
clear mpp


plot_oasis = 1;
trunc_oasis = 1;

count = 1; 

clear mpp
for i = 1:length(trials)
    
    cur_trial = trials(i);
    loc_id = loc_ids(i);
    this_seq{i} = data.trial_metadata(cur_trial).sequence;
    stims_per_trial(i) = length(this_seq{i});
    this_stim_key{i} = data.trial_metadata(cur_trial).stim_key;
    power_curve_num{i} = unique([this_seq{i}.target_power]);
    stim_starts{i} = [data.trial_metadata(cur_trial).sequence.start];
    
    if plot_oasis

        iter = iter_id(i);
        datafilename = [map_id '_z' num2str(loc_id) '_iter' num2str(iter) '.mat'];
        oasisfilename = [map_id '_z' num2str(loc_id) '_iter' num2str(iter) '_detect.mat'];

        load(datafilename)
        load(oasisfilename)
        oasis_data = reshape(event_process,size(traces'))'; 
    else
        oasis_data = zeros(length(this_seq{i}),1000);
    end
    
    for j = 1:length(this_seq{i})
        if i == 1 && j == 1
            full_seq(1) = this_seq{i}(j);
            if trunc_oasis
                mpp(1).times = ...
                    find(oasis_data(j,exp_data.params.time.min_time:exp_data.params.time.max_time),1)...
                     + exp_data.params.time.min_time - 1; 
            else
                mpp(1).times = ...
                    find(oasis_data(j,exp_data.params.time.min_time:exp_data.params.time.max_time))...
                     + exp_data.params.time.min_time - 1;
            end
        else
            full_seq(end+1) = this_seq{i}(j);
            if trunc_oasis
                mpp(end+1).times = ...
                    find(oasis_data(j,exp_data.params.time.min_time:exp_data.params.time.max_time),1)...
                     + exp_data.params.time.min_time - 1;
            else
                mpp(end+1).times = ...
                    find(oasis_data(j,exp_data.params.time.min_time:exp_data.params.time.max_time))...
                     + exp_data.params.time.min_time - 1;
            end
        end
        full_seq(end).precomputed_target_index = ...
            full_seq(end).precomputed_target_index + size(full_stim_key,1);
    end
    if size(this_stim_key{i},3) == 1
        this_stim_key{i} = cat(3,this_stim_key{i},nan([size(this_stim_key{i}) 2]));
    end
    full_stim_key = [full_stim_key; this_stim_key{i}];  
    
%     mpp(i)
end
power_curve_num = unique([power_curve_num{:}]);
maps = cell(length(power_curve_num),1);
this_seq = [this_seq{:}];
max_trial = length(this_seq);
% max_trial = 1200;
[traces_ch1,traces_ch2] = ...
    get_stim_stack(data,trials,...
        stims_per_trial,stim_starts,defaults.Fs,duration);
stim_inds = [full_seq.precomputed_target_index];
% on_cell_trials = isnan(full_stim_key(stim_inds,1,2));
on_cell_trials = ones(size(stim_inds))';
% power_curve_num = 150;
% traces = [];
stim_pow = [];
target_locs = [];
stim_inds = [];
deorder = [];
num_trials = 0;
spacing = 5;
% power_curve_num = power_curve_num(end-1:end);

mpp_pow = cell(length(power_curve_num),1);
% oasis_data = reshape(event_process,size(traces_ch1'))';
% mpp = struct();
% for i = 1:size(oasis_data,1)
%     mpp(i).times = find(oasis_data(i,45:300),1) + 44;
% end

for i= 1:length(full_seq)
    if full_seq(i).group == 3
        full_seq(i).target_power = 0;
    end
end

for i = 1%length(power_curve_num)
    
%     this_power = power_curve_num(i);
    this_power = 0;
    
%     this_seq = this_seq(1:max_trial);
    traces_pow{1} = traces_ch1(on_cell_trials' & [full_seq.target_power] == this_power,:);
%     traces = [traces; traces_pow{1}];
%     deorder = [deorder find(on_cell_trials' & [full_seq.target_power] == this_power)]; 
    traces_pow{2} = traces_ch2(on_cell_trials' & [full_seq.target_power] == this_power,:);
    this_seq_power = full_seq(on_cell_trials' & [full_seq.target_power] == this_power);
    mpp_pow{i}{1} = mpp(on_cell_trials' & [full_seq.target_power] == this_power);
    mpp_pow{i}{2} = mpp(on_cell_trials' & [full_seq.target_power] == this_power);
%     color_these_trials{i}{1} = group_colors([this_seq_power.group],:);
%     color_these_trials{i}{2} = group_colors([this_seq_power.group],:);
%     mpp_pow{i} = mpp(num_trials+(1:length(this_seq_power)));
%     mpp_pow{i}{1} = [];mpp_pow{i}{1} = [];
    num_trials = num_trials + length(this_seq_power);
    [maps_single{i}, mpp_maps{i}] = ...
...%         see_grid_multi(traces_pow,mpp_pow{i},this_seq_power,full_stim_key,spacing,1,0,0);
    see_grid_multi(traces_pow,mpp_pow{i},this_seq_power,full_stim_key,spacing,1,0,0);
%     title(['Power = ' num2str(power_curve_num(i)) ' mW'])
%     xlim(xlims); ylim(ylims);
%     get_mean_events_image(mpp_maps{i}, 2000, 1, 1);
%     title(['Event Counts, Power = ' num2str(power_curve_num(i)) ' mW'])
%     caxis([0 2])
end

%% gad2
xs = [10, 19, 22, 26, 25, 39, 56, 14, 17, 4, 32, 33]
ys = [53, 44, 54, 54, 58, 33, 43, 20, 22, 37, 42, 45]

%% emx site 1

xs = [5, 25, 30, 23]
ys = [52, 45, 30, 37]

%% emx site 2

xs = [24]
ys = [21]

%%
spots_use = 1:4;%[1, 4, 8, 10, 11, 12]
figure
count = 1;
for i = 1:length(spots_use)
    ind = spots_use(i)
    binx = xs(ind);
    biny = ys(ind);
    
    map_id = 3
%     figure
    subplot(length(spots_use),2,count)
    count = count + 1;
    plot_trace_stack(maps_single{map_id}{1}{binx,biny},0,'',[],[],[],1)
    ylim([-90 30])
    xlim([0 2000])
    subplot(length(spots_use),2,count)
    count = count + 1;
    if i == 4
        binx = 22; biny = 40;
    elseif i == 3
        binx = 29; biny = 30;
    end
    plot_trace_stack(maps_stp{map_id}{1}{binx,biny},0,'',[],[],[],1)
    ylim([-90 30])
    xlim([0 2000])
end



%%

these_cells = exp_data.cells_targets.cell_group_list{loc_id};
these_cell_locs = exp_data.cells_targets.cell_locations(these_cells,:);
