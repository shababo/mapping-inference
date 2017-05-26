function [pi_dense,inner_normalized_products] = get_weights_v2(...
    cell_params, shape_template,Z_dense)
% 
n_cell = size(cell_params.shape_gain,1);
pi_dense = zeros(n_cell,size(Z_dense,1));
for i_cell = 1:n_cell
    
    this_cell_shape = shape_template(cell_params.shape_gain(i_cell)).shape;
    % nucleus is assumed to be at the center
    center_idx=(size(this_cell_shape)+1)/2;
    % Relative location to the cell 
    relative_dist = Z_dense-ones(size(Z_dense,1),1)*cell_params.locations(i_cell,:);
    relative_idx = round(bsxfun(@plus,relative_dist,center_idx));
    
    
    for i_loc = 1:size(relative_idx,1)
        condition_shape = (relative_idx(i_loc,1) > 0 & relative_idx(i_loc,1) < size(this_cell_shape,1))& ...
            (relative_idx(i_loc,2) > 0 & relative_idx(i_loc,2) < size(this_cell_shape,2))& ...
            (relative_idx(i_loc,3) > 0 & relative_idx(i_loc,3) < size(this_cell_shape,3));
        if condition_shape
            
            pi_dense(i_cell,i_loc)=this_cell_shape(relative_idx(i_loc,1),relative_idx(i_loc,2),relative_idx(i_loc,3));
            if isnan(pi_dense(i_cell,i_loc))
                pi_dense(i_cell,i_loc)=0;
            end
        end
    end
end

%tend= toc;
%tend-tstart


weights_mat_full = pi_dense;
% Calculate the inner product of the induced probability
inner_products = weights_mat_full'*weights_mat_full;

self_products = diag(inner_products)*ones(1,size(inner_products,1));
inner_normalized_products = inner_products./self_products;


