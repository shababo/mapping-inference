function [this_neighbourhood] = undefined_to_connected(this_neighbourhood,group_profile)

group_ID='undefined';
i_cell_group_to_nhood= find(get_group_inds(this_neighbourhood,group_ID,this_neighbourhood.batch_ID));

switch group_profile.regroup_func_params.regroup_type
    case 'Quantiles'
        i_batch=this_neighbourhood.batch_ID;
        neurons=this_neighbourhood.neurons(i_cell_group_to_nhood);
        properties={'PR'};summary_stat={'lower_quantile','mean'};
        temp_output=grab_values_from_neurons(i_batch,neurons,properties,summary_stat);
        gamma_lower_quantile=temp_output.PR.lower_quantile;
        gamma_mean=temp_output.PR.mean;
        
        i_batch_prev=max(1,this_neighbourhood.batch_ID-1);
%         neurons=this_neighbourhood.neurons(i_cell_group_to_nhood);
        properties={'PR'};summary_stat={'mean'};
        temp_output=grab_values_from_neurons(i_batch_prev,neurons,properties,summary_stat);
        gamma_mean_previous=temp_output.PR.mean;
        
        max_changes_undefined= max(abs(gamma_mean_previous-gamma_mean)); %
        cell_list_undefined_to_connected = ...
           i_cell_group_to_nhood(gamma_lower_quantile > group_profile.regroup_func_params.disconnected_threshold);
    case 'NonzeroProb'
end

too_few_cells=(length(i_cell_group_to_nhood)-length(cell_list_undefined_to_connected)) < length(this_neighbourhood.neurons)*group_profile.regroup_func_params.singlespot_threshold;
too_tiny_change=max_changes_undefined<group_profile.regroup_func_params.change_threshold;

if too_few_cells || too_tiny_change
    cell_list_undefined_to_connected =i_cell_group_to_nhood;
end


if ~isempty(cell_list_undefined_to_connected)
    for i_cell = 1:length(cell_list_undefined_to_connected)
        this_cell=cell_list_undefined_to_connected(i_cell);
        this_neighbourhood.neurons(this_cell).group_ID{i_batch+1}='connected';
    end
end








