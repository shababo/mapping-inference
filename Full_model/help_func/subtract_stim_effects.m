function [prob_initial]=subtract_stim_effects(group_ID,this_cell,prob_initial,loc_selected, neurons)
    
    
    number_of_cells=length(neurons);
    pi_mat=[];
    for i_cell = 1:number_of_cells
        pi_mat=[pi_mat neurons(i_cell).stim_locations.(group_ID).effect(:,loc_selected(i_cell))];
        
    end
    pi_this=neurons(this_cell).stim_locations.(group_ID).effect(:,loc_selected(this_cell));
    inner_product_normalized= pi_mat'*pi_this/sqrt(pi_this'*pi_this);
    prob_initial=prob_initial-inner_product_normalized*prob_initial(this_cell);

end
