function [this_neighbourhood] = connected_to_alive(this_neighbourhood,group_profile)

group_ID='connected';
i_cell_group_to_nhood= find(get_group_inds(this_neighbourhood,group_ID));

switch group_profile.regroup_func_params.regroup_type
    case 'Quantiles'
        i_batch=this_neighbourhood.batch_ID;
        neurons=this_neighbourhood.neurons(i_cell_group_to_nhood);
        properties={'PR_params'};summary_stat={'lower_quantile','mean'};
        temp_output=grab_values_from_neurons(i_batch,neurons,properties,summary_stat);
        gamma_lower_quantile=temp_output.PR_params.lower_quantile;
        gamma_mean=temp_output.PR_params.mean;
        
        i_batch_prev=max(1,this_neighbourhood.batch_ID-1);
        %         neurons=this_neighbourhood.neurons(i_cell_group_to_nhood);
        properties={'PR_params'};summary_stat={'mean'};
        temp_output=grab_values_from_neurons(i_batch_prev,neurons,properties,summary_stat);
        gamma_mean_previous=temp_output.PR_params.mean;
        
        flags=   (abs(gamma_mean-gamma_mean_previous)<group_profile.regroup_func_params.change_threshold) & ...
            gamma_lower_quantile > group_profile.regroup_func_params.disconnected_threshold;
        cell_list_connected_to_alive = ...
                i_cell_group_to_nhood(flags);
        
    case 'NonzeroProb'
end


if ~isempty(cell_list_connected_to_alive )
    for i_cell = 1:length(cell_list_connected_to_alive )
        this_cell=cell_list_connected_to_alive(i_cell);
        this_neighbourhood.neurons(this_cell).group_ID{i_batch+1}='alive'; 
    end
end


end





