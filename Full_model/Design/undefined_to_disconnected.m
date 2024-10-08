function [this_neighbourhood] = undefined_to_disconnected(this_neighbourhood,group_profile)

group_ID='undefined';
i_cell_group_to_nhood= find(get_group_inds(this_neighbourhood,group_ID,this_neighbourhood.batch_ID));

switch group_profile.regroup_func_params.regroup_type
    case 'Quantiles'
        
        i_batch=this_neighbourhood.batch_ID;
        neurons=this_neighbourhood.neurons(i_cell_group_to_nhood);
        properties={'PR'};summary_stat={'upper_quantile'};
        temp_output=grab_values_from_neurons(i_batch,neurons,properties,summary_stat);

        gamma_upper_quantile=temp_output.PR.upper_quantile;
        cell_list_undefined_to_disconnected = ...
            i_cell_group_to_nhood(gamma_upper_quantile<group_profile.regroup_func_params.disconnected_threshold);
    case 'NonzeroProb'
end

if ~isempty(cell_list_undefined_to_disconnected)
    for i_cell = 1:length(cell_list_undefined_to_disconnected)
        this_cell=cell_list_undefined_to_disconnected(i_cell);
        this_neighbourhood.neurons(this_cell).group_ID{i_batch+1}='disconnected';
    end
end








