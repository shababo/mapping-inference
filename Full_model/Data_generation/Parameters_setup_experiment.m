%% Parameters for the experiments
%% Stimulus locations:

% Define a num_dense by num_dense grid
buffer=50;
x_dense = (0:(num_dense-1))*(2*buffer+max(Z(:,1))-min(Z(:,1)))/(num_dense-1) + min(Z(:,1))-buffer;
y_dense = (0:(num_dense-1))*(2*buffer+max(Z(:,2))-min(Z(:,2)))/(num_dense-1) + min(Z(:,2))-buffer;

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