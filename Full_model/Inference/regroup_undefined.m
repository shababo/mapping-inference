function [this_neighbourhood] = regroup_undefined(this_neighbourhood,to_groups,group_profile)

group_type_ID=group_profile.group_type_ID;
cells_this_group=index([this_neighbourhood.neurons(:).group_type_ID]==group_type_ID);
number_cells_this_group=length(cells_this_group);
number_cells_all= length(this_neighbourhood.neurons);



undefined_profile.regroup_function=@regroup_undefined;
undefined_profile.regroup_func_params=struct;
undefined_profile.regroup_func_params.connected_threshold=0.5;
undefined_profile.regroup_func_params.disconnected_threshold=0.2;
undefined_profile.regroup_func_params.quantile_prob=[0.05 0.95];
undefined_profile.regroup_func_params.regroup_type='Quantiles'; % Quantiles or NonzeroProb
undefined_profile.regroup_func_params.singlespot_threshold=0.2;% certain proportion of cells 
undefined_profile.regroup_func_params.change_threshold =0.05;

switch regroup_type
    case 'Quantiles'
        gamma_upper_quantile=extract_this_entry(this_neighbourhood.neurons(cells_this_group));
        gamma_lower_quantile=extract_this_entry(this_neighbourhood.neurons(cells_this_group));
        gamma_mean=extract_this_entry(this_neighbourhood.neurons(cells_this_group));
        gamma_mean_previous=extract_this_entry(this_neighbourhood.neurons(cells_this_group)); 
        max_changes_undefined= max(abs(gamma_mean_previous-gamma_mean)); %
        undefined_to_disconnected = ...
            intersect(find([parameter_path(iter+1,:).gamma_upper_quantile]<group_profile.regroup_func_params.disconnected_threshold),...
            find( undefined_cells_current));
        undefined_to_connected = ...
            intersect(find([parameter_path(iter+1,:).gamma_mean]>group_profile.regroup_func_params.connected_threshold
        ),find( undefined_cells_current));
        
        
    case 'NonzeroProb'
end
if (sum(output.dead_cells)+sum(output.alive_cells))  == length(output.alive_cells)
    output.connected_cells(find(output.alive_cells))=1;

    
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





