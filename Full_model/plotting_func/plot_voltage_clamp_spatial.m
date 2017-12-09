function plot_voltage_clamp_spatial(experiment_setup,neighbourhood_ID,batch_IDs,varargin)

if ~isempty(varargin) && ~isempty(varargin{1})
    
    dirpath = varargin{1};
 
else
    
    dirpath = experiment_setup.analysis_root;
end
    
trials = struct();
group_names = experiment_setup.group_names;

% load batches and collect trials
for i = 1:length(batch_IDs)
    
    batchsavepath = [dirpath experiment_setup.exp_id ...
                    '_n' num2str(neighbourhood_ID)...
                    '_b' num2str(batch_IDs(i)) '_complete.mat'];
                
    load(batchsavepath)
    for j = 1:length(group_names)

        trials = [trials experiment_query.(group_names{j}).trials];
        
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