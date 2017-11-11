function [pi_stim] = get_weights(...
    cell_params, template_cell,stim_locations)

% NEED TO ALLOW FOR CELL-SPECIFIC EFFECT

number_of_cells = length(cell_params);
pi_stim = zeros(number_of_cells,size(stim_locations,1));
for i_cell = 1:number_of_cells 
    
    this_cell_shape = template_cell.cell_shape;
    % nucleus is assumed to be at the center
    center_idx=(size(this_cell_shape)+1)/2;
    % Relative location to the cell
    relative_dist = stim_locations-ones(size(stim_locations,1),1)*cell_params(i_cell).location;
    relative_idx = round(bsxfun(@plus,relative_dist,center_idx));
    for i_loc = 1:size(relative_idx,1)
        condition_shape = (relative_idx(i_loc,1) > 0 & relative_idx(i_loc,1) < size(this_cell_shape,1))& ...
            (relative_idx(i_loc,2) > 0 & relative_idx(i_loc,2) < size(this_cell_shape,2))& ...
            (relative_idx(i_loc,3) > 0 & relative_idx(i_loc,3) < size(this_cell_shape,3));
        if condition_shape
            pi_stim(i_cell,i_loc)=this_cell_shape(relative_idx(i_loc,1),relative_idx(i_loc,2),relative_idx(i_loc,3));
            if isnan(pi_stim(i_cell,i_loc))
                pi_stim(i_cell,i_loc)=0;
            end
        end
    end
end

% Calculate the inner product of the induced probability
% inner_products = pi_stim'*pi_stim;
% self_products = diag(inner_products)*ones(1,size(inner_products,1));
% inner_normalized_products = inner_products./self_products;


