function [pilot_data] = map_result_to_neuron(path,varargin)

if ~isempty(varargin)
    oldpath =varargin{1};
end
%%
axis_list= {'x' 'y' 'z' 'xy'};
pilot_data = struct;


for i_axis = 1:4
    target_axis=axis_list{i_axis};
    ax=target_axis;
      clear('neurons')
          
    switch target_axis
        case 'z'
            load(path.(ax));
            n_cell = length(result_z);
            neurons(n_cell)=struct;
            for i_cell = 1:n_cell
                neurons(i_cell).stim_grid= reshape(result_z(i_cell).current_z_pos, 1,[]); % row vector
                %       neurons(i_cell).max_current=result_xy(i_cell).(ax).current';
                neurons(i_cell).power=50*ones(1,length(neurons(i_cell).stim_grid));
                neurons(i_cell).raw_current=result_z(i_cell).max_curr';
            end
            
        case 'x'
            
            load(oldpath);
            %             result_xy=result_xy;
            n_cell = length(result_xy);
            % Say we are interested in x-shape
            % Count the number of trials with the same y-coordinate
            % The y-coordinate with the largest number of trials is then y0
            for i_cell = 1:length(result_xy)
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
                %                 result_xy(i_cell).(ax).spike_times = result_xy(i_cell).spike_times(this_axis_trials);
            end
            if isfield(path,ax)
                load(path.(ax));
                new_cell = length(result_x);
                neurons(n_cell+new_cell)=struct;
            else
                neurons(n_cell)=struct;
            end
            for i_cell = 1:n_cell
                neurons(i_cell).stim_grid= result_xy(i_cell).(target_axis).stim_grid';
                %       neurons(i_cell).max_current=result_xy(i_cell).(ax).current';
                neurons(i_cell).power=result_xy(i_cell).(target_axis).power';
                neurons(i_cell).raw_current=result_xy(i_cell).(target_axis).current';
                %             neurons(i_cell).spike_times=result_xy(i_cell).(target_axis).spike_times';
                %         neurons(i_cell).scaled_current= neurons(i_cell).scaled_current/max( neurons(i_cell).scaled_current);
            end
            if isfield(path,ax)
                
                for i_cell = 1:new_cell
                    neurons(i_cell+n_cell).stim_grid= reshape(result_x(i_cell).current_x_pos, 1,[]); % column vector
                    neurons(i_cell+n_cell).power=50*ones(1,length(neurons(i_cell+n_cell).stim_grid));
                    neurons(i_cell+n_cell).raw_current=result_x(i_cell).max_curr';
                end
            end
        case 'y'
            load(oldpath);
            %             result_xy=result_xy_nuc_detect;
            n_cell = length(result_xy);
            % Say we are interested in x-shape
            % Count the number of trials with the same y-coordinate
            % The y-coordinate with the largest number of trials is then y0
            for i_cell = 1:length(result_xy)
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
                %                 result_xy(i_cell).(ax).spike_times = result_xy(i_cell).spike_times(this_axis_trials);
            end
            result_xy(4)=result_xy(2);
            result_xy(7) = result_xy(2);
            
            if isfield(path,ax)
                load(path.(ax));
                new_cell = length(result_x);
                neurons(n_cell+new_cell)=struct;
            else
                neurons(n_cell)=struct;
            end
            for i_cell = 1:n_cell
                neurons(i_cell).stim_grid= result_xy(i_cell).(target_axis).stim_grid';
                %       neurons(i_cell).max_current=result_xy(i_cell).(ax).current';
                neurons(i_cell).power=result_xy(i_cell).(target_axis).power';
                neurons(i_cell).raw_current=result_xy(i_cell).(target_axis).current';
                %             neurons(i_cell).spike_times=result_xy(i_cell).(target_axis).spike_times';
                %         neurons(i_cell).scaled_current= neurons(i_cell).scaled_current/max( neurons(i_cell).scaled_current);
            end
            if isfield(path,ax)
                
                for i_cell = 1:new_cell
                    neurons(i_cell+n_cell).stim_grid= reshape(result_x(i_cell).current_x_pos, 1,[]); % column vector
                    neurons(i_cell+n_cell).power=50*ones(1,length(neurons(i_cell+n_cell).stim_grid));
                    neurons(i_cell+n_cell).raw_current=result_x(i_cell).max_curr';
                end
            end
        case 'xy'
            if isfield(path,ax)
                load(path.(ax));
                n_cell = length(result_shape);
                neurons(n_cell)=struct;
                for i_cell = 1:n_cell
                    % Only take the spots at 200 microns
                    center_indices=find(result_shape(i_cell).current_targ_pos(:,3)==200);
                    neurons(i_cell).stim_grid= result_shape(i_cell).current_targ_pos(center_indices,:)';
                    %       neurons(i_cell).max_current=result_xy(i_cell).(ax).current';
                    neurons(i_cell).power=50*ones(1,size(neurons(i_cell).stim_grid,2));
                    neurons(i_cell).raw_current=result_shape(i_cell).max_curr(center_indices)';
                end
            end
            
    for i_cell = 1:n_cell
        figure(i_cell)
        hist(result_shape(i_cell).max_curr)
    end
    end
    if exist('neurons')
    pilot_data.(target_axis)=struct;
    pilot_data.(target_axis).neurons=neurons;
    end
    
end
end
