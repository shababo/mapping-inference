function [pi_dense,inner_normalized_products] = get_weights(...
    A,locations,shape_mat,threshold,Z_dense)
%% Pre-calculate the spatial marks as the stimulus convoluted with the cell shapes
% Calculate the distance to the center of the stimulus:
% Sum up the total effects for each cell 
% Calculate the probability of firing for ALL neurons
% We will use it in simulating the real spikes
%tic;
%tstart = toc;
n_cell = size(locations,1);
pi_dense = zeros(n_cell,size(Z_dense,1));
for i_cell = 1:n_cell
    center_dist = max(exp(-0.5*pdist2(locations(i_cell,:),Z_dense,'mahalanobis',A).^2));
    % normalizing constant 
    diff_mat_cell = pdist2(shape_mat{i_cell},locations(i_cell,:),'mahalanobis',A);
    pi_mat = exp(-0.5*(diff_mat_cell.^2));
    pi_norm= sum(pi_mat,1);
    if center_dist > threshold
        diff_mat_cell = pdist2(shape_mat{i_cell},Z_dense,'mahalanobis',A);
        pi_dense_cell = exp(-0.5*(diff_mat_cell.^2));
        pi_dense(i_cell,:) = sum(pi_dense_cell,1)/pi_norm;
    end
end

%tend= toc;
%tend-tstart


weights_mat_full = pi_dense;
% Calculate the inner product of the induced probability
inner_products = weights_mat_full'*weights_mat_full;

self_products = diag(inner_products)*ones(1,size(inner_products,1));
inner_normalized_products = inner_products./self_products;


