function [neurons] = map_result_to_neuron(target_axis, path)

switch target_axis
    case 'z'
        load(path);
        n_cell = length(result_z);
        clear('neurons')
        neurons(n_cell)=struct;
        for i_cell = 1:n_cell
            neurons(i_cell).stim_grid= reshape(result_z(i_cell).current_z_pos, 1,[]); % row vector
            %       neurons(i_cell).max_current=result_xy(i_cell).(ax).current';
            neurons(i_cell).power=50*ones(1,length(neurons(i_cell).stim_grid));
            neurons(i_cell).raw_current=result_z(i_cell).max_curr';
        end
    case 'x'
        load(path);        
        n_cell = length(result_x);
        clear('neurons')
        neurons(n_cell)=struct;
        for i_cell = 1:n_cell
            neurons(i_cell).stim_grid= reshape(result_x(i_cell).current_x_pos, 1,[]); % column vector
            neurons(i_cell).power=50*ones(1,length(neurons(i_cell).stim_grid));
            neurons(i_cell).raw_current=result_x(i_cell).max_curr';
        end
end

