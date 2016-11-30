
%% Option II: fit a full model but use the online update 
% This might be slow
% Good news is that it does not depend on Stage I too much 
% To improve the effiiency, we need to reduce the density of the grid 
%% Design parameters:

%% For random stimuli
%% Find the neurons that locate in the region:
% 
% all_neuron_locations record locations of the neurons and their layer info
all_neuron_locations =zeros(0,4);
for i = 1:num_layers
    all_neuron_locations = [all_neuron_locations; [neuron_locations{i}(:,1:3) i*ones(size(neuron_locations{i},1),1) ]];
end

epsilon = 0.01;
x_low = postsyn_position(1) - region_width/2-epsilon;
x_upp = postsyn_position(1) + region_width/2;
y_low = postsyn_position(2) - region_height/2-epsilon;
y_upp = postsyn_position(2) + region_height/2;

% Identify neurons within this location:
neuron_in_region = zeros(size(all_neuron_locations,1),1); 
for i = 1:size(all_neuron_locations,1)
    if all_neuron_locations(i,1) > x_low & all_neuron_locations(i,1) < x_upp
        if all_neuron_locations(i,2) > y_low & all_neuron_locations(i,2) < y_upp
            neuron_in_region(i)=1;
        end
    end
end

Z = zeros(sum(neuron_in_region),3);
count = 1;
for i = 1:size(all_neuron_locations,1)
    if neuron_in_region(i) > 0
        Z(count,:) = [all_neuron_locations(i,1:2) postsyn_position(3)];
        count = count + 1;
    end 
end


K_z = size(Z,1);

B = diag(all_locations*inv(A)*all_locations');
diff_mat=B*ones(size(Z,1),1)'+(ones(size(all_locations,1),1)*diag(Z*inv(A)*Z')')-2*all_locations*inv(A)*Z';
pi_Z_all = exp(-0.5*(diff_mat));

% Calculate the probability of firing for the LOCAL neurons (i.e., those
% that are within the 2-D plane we consider)
% We will use this in estimating the expectation, and in fitting the model
B = diag(Z*inv(A)*Z');
diff_mat=B*ones(size(Z,1),1)'+(ones(size(Z,1),1)*diag(Z*inv(A)*Z')')-2*Z*inv(A)*Z';
pi_Z_grid = exp(-0.5*(diff_mat));

%% Define the combination of trials 
% We want to stimulate neurons that are far apart from each other

% To do this, we split the cells by the quantiles of their x and y
% coordinates, then alternating between x and y.

% Find the quantiles:


[~, ~, x_ranks] = unique(Z(:,1));
x_freq = x_ranks / size(Z,1);

x_index = ceil(x_freq / (1/num_sources));

x_group = cell(num_sources,1);
for i = 1:size(Z,1)
    x_group{x_index(i)} = [x_group{x_index(i)} i];
end
max_x_group = 0;
for i = 1:num_sources;
    max_x_group = max(max_x_group, size(x_group{i},2));
end

% In the simulation, select one site from each group 


%% For entropy estimation
% The new stimuli are chosen to stay away from the true neurons
% We stimulate a dense grid instead
%x_dense = zeros(num_dense,1);
%y_dense = zeros(num_dense,1);
x_dense = (0:(num_dense-1))*(max(Z(:,1))-min(Z(:,1)))/(num_dense-1) + min(Z(:,1));
y_dense = (0:(num_dense-1))*(max(Z(:,2))-min(Z(:,2)))/(num_dense-1) + min(Z(:,2));
Z_dense = zeros(num_dense^2,3);
for i = 1:num_dense
    for l = 1:num_dense
        Z_dense((i-1)*num_dense + l,:) = [x_dense(i) y_dense(l) postsyn_position(3)];
    end
end

% Calculate the probability of firing for ALL neurons
% We will use it in simulating the real spikes
B = diag(all_locations*inv(A)*all_locations');
diff_mat=B*ones(size(Z_dense,1),1)'+(ones(size(all_locations,1),1)*diag(Z_dense*inv(A)*Z_dense')')-2*all_locations*inv(A)*Z_dense';
pi_dense_all = exp(-0.5*(diff_mat));

% Calculate the probability of firing for the LOCAL neurons (i.e., those
% that are within the 2-D plane we consider)
% We will use this in estimating the expectation, and in fitting the model
B = diag(Z*inv(A)*Z');
diff_mat=B*ones(size(Z_dense,1),1)'+(ones(size(Z,1),1)*diag(Z_dense*inv(A)*Z_dense')')-2*Z*inv(A)*Z_dense';
pi_dense_grid = exp(-0.5*(diff_mat));
%% Calculate and save the innner products of weight matrix w
% The responsive cells might not cover the full space
% Hence, we can further save some computation by considering only the
% locations that are relevant (with total firing probability > 0.1)
weights_mat_full = pi_dense_grid;

% Calculate the inner product of the induced probability
inner_products = weights_mat_full'*weights_mat_full;

self_products = diag(inner_products)*ones(1,size(inner_products,1));
inner_normalized_products = inner_products./self_products;
%
% a= weights_mat_full(:,1);
% b= weights_mat_full(:,2);
%
% k= a'*b/b'*b
% a'*a /b'*b + k^2  - 2k a'*b/b'*b
% (a- k*b)'*(a- k*b)
