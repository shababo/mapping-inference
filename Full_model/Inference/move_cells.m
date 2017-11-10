function [output] = move_cells(...
    connected_threshold,disconnected_threshold,single_spot_threshold,multispot_change_threshold,...
    dead_cells_current,disconnected_cells_current,undefined_cells_current,potentially_connected_cells_current,alive_cells_current,...
    variational_params_path, parameter_path,assignment_type)

output=struct;
output.undefined_cells=undefined_cells_current;
output.connected_cells=potentially_connected_cells_current;
output.disconnected_cells=disconnected_cells_current;
output.dead_cells=dead_cells_current;
output.alive_cells=alive_cells_current;

if (sum(output.dead_cells)+sum(output.alive_cells))  == length(output.alive_cells)
    output.connected_cells(find(output.alive_cells))=1;
else
    
    %------------------------------------------------------%
    % Moving the cell between groups
    % mean_gamma_undefined & undefined_cells{iter} % A
    % mean_gamma_disconnected & potentially_disconnected_cells{iter} %B
    % mean_gamma_connected & potentially_connected_cells{iter} %C
    undefined_cell_list=find(undefined_cells_current);
    iter=size(parameter_path,1)-1;
    switch assignment_type
        case 1 % gamma quantiles
            undefined_to_disconnected = ...
                intersect(find([parameter_path(iter+1,:).gamma_upper_quantile]<disconnected_threshold),...
                find( undefined_cells_current));
            connected_to_disconnected = ...
                intersect(find([parameter_path(iter+1,:).gamma_upper_quantile]<disconnected_threshold),...
                find(potentially_connected_cells_current));
            undefined_to_connected = ...
                intersect(find([parameter_path(iter+1,:).gamma_mean]>connected_threshold),...
                find( undefined_cells_current));
            connected_to_alive = ...
                intersect(find([parameter_path(iter+1,:).gamma_lower_quantile]> disconnected_threshold),...
                find(potentially_connected_cells_current));
            max_changes_undefined=...
                max(abs([parameter_path(iter+1,undefined_cell_list).gamma_mean]-[parameter_path(iter,undefined_cell_list).gamma_mean])); %
            %+...
            %    max(abs([parameter_path(iter+1,undefined_cell_list).gain_mean]-[parameter_path(iter,undefined_cell_list).gain_mean])./...
            %   [parameter_path(iter,undefined_cell_list).gain_mean]);
        case 2 % gamma and spike probability
            nonzero_prob = [parameter_path(iter+1,:).nonzero_prob];
            
            undefined_to_disconnected = ...
                intersect(find(nonzero_prob <disconnected_threshold),...
                find( undefined_cells_current));
            connected_to_disconnected = ...
                intersect(find(nonzero_prob <disconnected_threshold),...
                find(potentially_connected_cells_current));
            undefined_to_connected = ...
                intersect(find(nonzero_prob >connected_threshold),...
                find( undefined_cells_current));
            connected_to_alive = ...
                intersect(find(nonzero_prob > connected_threshold),...
                find(potentially_connected_cells_current));
            max_changes_undefined=...
                max(abs(nonzero_prob -[parameter_path(iter,:).nonzero_prob]));
            
            %+...
            %   max(abs([parameter_path(iter+1,undefined_cell_list).gain_mean]-[parameter_path(iter,undefined_cell_list).gain_mean])./...
            %  [parameter_path(iter,undefined_cell_list).gain_mean]);
    end
    %     if synchrony_test
    %        disconnected_to_dead=find(disconnected_indicators & disconnected_cells_current);
    %        disconnected_to_connected=find( (~disconnected_indicators) & disconnected_cells_current);
    %
    %     else
    %        disconnected_to_dead=find(disconnected_cells_current);
    %        disconnected_to_connected=([]);
    %     end
    
    % Declare a cell to be disconnected if its posterior variances of gains
    % are large, or its posterior spike probability is also
    % Update the cell lists:
    output=struct;
    output.undefined_cells=undefined_cells_current;
    output.undefined_cells(undefined_to_disconnected)=0;
    output.undefined_cells(undefined_to_connected)=0;
    
    % after moving cells_out
    if sum(output.undefined_cells)<single_spot_threshold ||  max_changes_undefined<multispot_change_threshold
        undefined_to_connected = [undefined_to_connected;  find(output.undefined_cells)];
    end
    output.undefined_cells(undefined_to_connected)=0;
    
    output.connected_cells=potentially_connected_cells_current;
    output.connected_cells(connected_to_disconnected)=0;
    output.connected_cells(undefined_to_connected)=1;
    output.connected_cells(connected_to_alive)=0; % keep collecting data on these cells
    %     output.connected_cells(disconnected_to_connected)=1; % keep collecting data on these cells
    
    %     potentially_connected_cells{iter+1}(connected_to_alive)=0;
    
    output.disconnected_cells=disconnected_cells_current;
    output.disconnected_cells(connected_to_disconnected)=1;
    output.disconnected_cells(undefined_to_disconnected)=1;
    %     output.disconnected_cells(disconnected_to_dead)=0;
    %     output.disconnected_cells(disconnected_to_connected)=0;
    
    output.dead_cells=dead_cells_current;
    %     output.dead_cells(disconnected_to_dead)=1;
    
    output.alive_cells=alive_cells_current;
    output.alive_cells(connected_to_alive)=1;
    
end

end





