function [this_neighbourhood] = connected_to_alive(this_neighbourhood,group_profile)

group_ID='connected';
cells_this_group= find(get_group_inds(this_neighbourhood,group_ID));

switch group_profile.regroup_func_params.regroup_type
    case 'Quantiles'
           i_batch=this_neighbourhood.batch_ID;
        neurons=this_neighbourhood.neurons(cells_this_group);
        properties={'PR_params'};summary_stat={'lower_quantile','mean'};
        temp_output=grab_values_from_neurons(i_batch,neurons,properties,summary_stat);
        gamma_lower_quantile=temp_output.PR_params.lower_quantile;
        gamma_mean=temp_output.PR_params.mean;
        
        cell_list_connected_to_alive = ...
                cells_this_group(gamma_lower_quantile> group_profile.regroup_func_params.disconnected_threshold);
        
    case 'NonzeroProb'
end




if ~isempty(cell_list_connected_to_alive )
    for i_cell = 1:length(cell_list_connected_to_alive )
        this_cell=cell_list_connected_to_alive (i_cell);
        this_neighbourhood.neurons(this_cell).group_ID='alive';
    end
end


end





