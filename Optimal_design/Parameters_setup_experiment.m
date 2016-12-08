%% Parameters for the experiments
%% Stimulus locations:

% Define a num_dense by num_dense grid
x_dense = (0:(num_dense-1))*(max(Z(:,1))-min(Z(:,1)))/(num_dense-1) + min(Z(:,1));
y_dense = (0:(num_dense-1))*(max(Z(:,2))-min(Z(:,2)))/(num_dense-1) + min(Z(:,2));
Z_dense = zeros(num_dense^2,3);
for i = 1:num_dense
    for l = 1:num_dense
        Z_dense((i-1)*num_dense + l,:) = [x_dense(i) y_dense(l) postsyn_position(3)];
    end
end

% Merge the grid with the centers of neurons 
Z_dense = [Z_dense; Z];

%% Pre-calculate the induced firing probability 

% Calculate the probability of firing for ALL neurons
% We will use it in simulating the real spikes
diff_mat = pdist2(all_locations,Z_dense,'mahalanobis',A);
pi_dense_all = exp(-0.5*(diff_mat.^2));

% Calculate the probability of firing for the LOCAL neurons (i.e., those
% that are within the 2-D plane we consider)
% We will use this in estimating the expectation, and in fitting the model
diff_mat = pdist2(Z,Z_dense,'mahalanobis',A);
pi_dense_local = exp(-0.5*(diff_mat.^2));

%% Calculate and save the innner products of weight matrix w

weights_mat_full = pi_dense_local;
% Calculate the inner product of the induced probability
inner_products = weights_mat_full'*weights_mat_full;

self_products = diag(inner_products)*ones(1,size(inner_products,1));
inner_normalized_products = inner_products./self_products;
%
% a= weights_mat_full(:,1);
% b= weights_mat_full(:,2);
% 
% k= (a'*b)/(b'*b)
% (a'*a)/(b'*b) 
% k^2  
% 2*k*(a'*b)/(b'*b)
% 
% (a- k*b)'*(a- k*b)
% (a- k*b)'*b

%% Put cells into different groups (the initialization/orthogonal design)
% We want to stimulate neurons that are far apart from each other
nquantile = num_sources*2;
[~, ~, x_ranks] = unique(Z(:,1));
x_freq = x_ranks / size(Z,1);
x_index = ceil(x_freq /( max(x_freq)/nquantile) );
x_group = cell(nquantile,1);
for i = 1:size(Z,1)
    x_group{x_index(i)} = [x_group{x_index(i)} i];
end

max_x_group = 0;
for i = 1:nquantile
    %size(x_group{i},2)
    max_x_group = max(max_x_group, size(x_group{i},2));
end
for i = 1:nquantile
    if size(x_group{i},2) < max_x_group
        x_group{i} = [x_group{i} randsample(x_group{i}, max_x_group-size(x_group{i},2))];
    end
end

diff_mat = pdist2(all_locations,Z,'mahalanobis',A);
pi_Z_all = exp(-0.5*(diff_mat.^2));

diff_mat = pdist2(Z,Z,'mahalanobis',A);
pi_Z_local = exp(-0.5*(diff_mat.^2));
