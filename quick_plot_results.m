
% close all
trials = 3:8;
loc_id = 2;
% trials = 4;

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

%%

group_colors = [0 0 1;
                1 .8 .8;
                .8 1 .8;
                1 0 0;
                0 1 0];
            
group_names = {'undefined_cells','potentially_disconnected_cells',...
                'potentially_connected_cells','dead_cells','alive_cells'};
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
%%

group_colors = [0 0 1;
                1 .8 .8;
                .8 1 .8;
                1 0 0;
                0 1 0];
            
group_names = {'undefined_cells','potentially_disconnected_cells',...
                'potentially_connected_cells','dead_cells','alive_cells'};
iter = length(exp_data.design.(group_names{1}){loc_id});
figure
for i = 1:25
    subplot(5,5,i)
     iter = i
    for i = 1:length(group_names) 
    
        
       
    these_cells = logical(exp_data.design.(group_names{i}){loc_id}{iter});
    scatter(cell_group_locs(these_cells,2),-cell_group_locs(these_cells,1),...
        exp_data.design.gamma_path{loc_id}(these_cells,iter)*20,...
        repmat(group_colors(i,:),sum(these_cells),1),...
        'filled');
    hold on
    end
end
%%

this_seq = cell(length(trials),1);
this_stim_key = cell(length(trials),1);
power_curve_num = cell(length(trials),1);
stim_starts = cell(length(trials),1);

full_stim_key = [];
clear full_seq
% full_seq = struct();
count = 1;
for i = 1:length(trials)
    cur_trial = trials(i);
    this_seq{i} = data.trial_metadata(cur_trial).sequence;
    stims_per_trial(i) = length(this_seq{i});
    this_stim_key{i} = data.trial_metadata(cur_trial).stim_key;
    power_curve_num{i} = unique([this_seq{i}.target_power]);
    stim_starts{i} = [data.trial_metadata(cur_trial).sequence.start];
    for j = 1:length(this_seq{i})
        if i == 1 && j == 1
            full_seq(1) = this_seq{i}(j);
        else
            full_seq(end+1) = this_seq{i}(j);
        end
        full_seq(end).precomputed_target_index = ...
            full_seq(end).precomputed_target_index + size(full_stim_key,1);
    end
    if size(this_stim_key{i},3) == 1
        this_stim_key{i} = cat(3,this_stim_key{i},nan([size(this_stim_key{i}) 2]));
    end
    full_stim_key = [full_stim_key; this_stim_key{i}];
    mpp(i)
end
power_curve_num = unique([power_curve_num{:}]);
maps = cell(length(power_curve_num),1);
this_seq = [this_seq{:}];
max_trial = length(this_seq);
% max_trial = 1200;
[traces_ch1,traces_ch2] = ...
    get_stim_stack(data,trials,...
        stims_per_trial,stim_starts,defaults.Fs);
stim_inds = [full_seq.precomputed_target_index];
% on_cell_trials = isnan(full_stim_key(stim_inds,1,2));
on_cell_trials = ones(size(stim_inds))';
% power_curve_num = 150;
traces = [];
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

for i = 1:length(power_curve_num)
    
    
%     this_seq = this_seq(1:max_trial);
    traces_pow{1} = traces_ch1(on_cell_trials' & [full_seq.target_power] == power_curve_num(i),:);
%     traces = [traces; traces_pow{1}];
%     deorder = [deorder find(on_cell_trials' & [full_seq.target_power] == power_curve_num(i))]; 
    traces_pow{2} = traces_ch2(on_cell_trials' & [full_seq.target_power] == power_curve_num(i),:);
    this_seq_power = full_seq(on_cell_trials' & [full_seq.target_power] == power_curve_num(i));
%     mpp_pow{i} = mpp(on_cell_trials' & [full_seq.target_power] == power_curve_num(i));
%     mpp_pow{i} = mpp(num_trials+(1:length(this_seq_power)));
    mpp_pow{i} = [];
    num_trials = num_trials + length(this_seq_power);
    [maps{i}, mpp_maps{i}] = see_grid_multi(traces_pow,mpp_pow{i},this_seq_power,full_stim_key,spacing,1,1);
%     title(['Power = ' num2str(power_curve_num(i)) ' mW'])
%     xlim(xlims); ylim(ylims);
%     get_mean_events_image(mpp_maps{i}, 2000, 1, 1);
%     title(['Event Counts, Power = ' num2str(power_curve_num(i)) ' mW'])
%     caxis([0 2])
end

%%

these_cells = exp_data.cells_targets.cell_group_list{loc_id};
these_cell_locs = exp_data.cells_targets.cell_locations(these_cells,:);