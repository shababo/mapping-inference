function [neurons, trials] = map_result_to_mapping_struct(target_axis, z_path, xy_path)
%%
load(z_path);
result_z(1)=result_z(2);
load(xy_path);
result_xy=result_xy_correct_power;

axis_list={'x' 'y' 'z'};
good_cell_list= [1:3 5:length(result_xy)];

for i_cell = good_cell_list
    
    for i_axis = 1:2
        ax=axis_list{i_axis};
        
        % Say we are interested in x-shape
        % Count the number of trials with the same y-coordinate
        % The y-coordinate with the largest number of trials is then y0
        other_idx = 3-i_axis; % 2 or 1
        [unique_grid,ia,ic]=unique(abs(result_xy(i_cell).current_targ_pos(:,other_idx)));
        a_counts = accumarray(ic,1);
        [~,most_freq_loc]=max(a_counts);
        loc=unique_grid(most_freq_loc); % 
        this_axis_trials = find(ic==most_freq_loc);
        
        % Now construct a structure for these trials:
        result_xy(i_cell).(ax)=struct;
        result_xy(i_cell).(ax).trials = this_axis_trials;
        result_xy(i_cell).(ax).stim_grid = result_xy(i_cell).current_targ_pos(this_axis_trials,i_axis);
        result_xy(i_cell).(ax).current = result_xy(i_cell).max_curr(this_axis_trials);
        result_xy(i_cell).(ax).power = result_xy(i_cell).curr_targ_power(this_axis_trials);
        result_xy(i_cell).(ax).spike_times = result_xy(i_cell).spike_times(this_axis_trials);
        
        
    end
    
end
% replace two bad cells with cell 2;
result_xy(4)=result_xy(2); 
% result_xy(7) = result_xy(2); 

i_ax=find(strcmp(target_axis,axis_list));

if ~strcmp(target_axis,'z')
    n_cell = length(result_xy);
    clear('neurons_tmp')
    neurons_tmp(n_cell)=struct;
    for i_cell = 1:n_cell
        neurons_tmp(i_cell).stim_grid= result_xy(i_cell).(target_axis).stim_grid';
        %       neurons(i_cell).max_current=result_xy(i_cell).(ax).current';
        neurons_tmp(i_cell).power=result_xy(i_cell).(target_axis).power';
        
        neurons_tmp(i_cell).raw_current=result_xy(i_cell).(target_axis).current';
        neurons_tmp(i_cell).spike_times=result_xy(i_cell).(target_axis).spike_times';
        
        %         neurons(i_cell).scaled_current= neurons(i_cell).scaled_current/max( neurons(i_cell).scaled_current);
    end
else
    n_cell = length(result_z);
    clear('neurons')
    neurons_tmp(n_cell)=struct;
    
    for i_cell = 1:n_cell
        neurons_tmp(i_cell).stim_grid=result_z(i_cell).current_z_pos;
        neurons_tmp(i_cell).raw_current=result_z(i_cell).max_curr';
        neurons_tmp(i_cell).power=50*ones(1,length(neurons_tmp(i_cell).stim_grid));
        neurons_tmp(i_cell).spike_times=result_z(i_cell).spike_times';
        % power is the same across trials
    end
    
end

%%
tmp_params=initialize_prior(['C:/Users/Shizhe/Documents/Mapping_data/Data/' 'Summary180627.mat']);
tmp_params=rmfield(tmp_params,{'GP_params' 'boundary_params','shift_y','shift_z'});
   

%% Generate the neurons & trials structures
total_cell=n_cell;
total_trials = length([neurons_tmp(:).power]);
trials(total_trials)=struct;
neurons(n_cell)=struct;
max_gap = 100;quantile_prob=[0.05 0.5 0.95];
i_trial = 0;
for i_cell = 1:total_cell
    neurons(i_cell).params(1) =tmp_params;
    neurons(i_cell).location=i_cell*max_gap;
     neurons(i_cell).posterior_stat(1)=...
            calculate_posterior(neurons(i_cell).params(1),...
            quantile_prob);
    for i = 1:length(neurons_tmp(i_cell).power)
        i_trial = i_trial +1;
        trials(i_trial).locations = neurons_tmp(i_cell).stim_grid(i)+neurons(i_cell).location;
        trials(i_trial).power_levels =  neurons_tmp(i_cell).power(i);
        trials(i_trial).event_times =  neurons_tmp(i_cell).spike_times(i);
        
    end
        
end
%%