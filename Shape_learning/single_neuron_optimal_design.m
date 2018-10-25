function [oracle_locations,optimal_locations ]=single_neuron_optimal_design(one_neuron,design_params, params,prior_params, prior_info)
% 

% sigma_epsilon =  0.1;
x_grid_sample =-30:2:30;
y_grid_sample = -30:2:30;
  S=100;
n_grid =51;
this_z =0;
n_locations = design_params.nlocs_per_neuron;
gap= design_params.gap;
%% True shape:
  % Calculate the shape values at the mesh grid (similar to visualize_3D_GP

this_z_indices=find(params.mesh_grid(:,3)==this_z);
true_shape_values=one_neuron.truth.shape(this_z_indices);
xy_locations=params.mesh_grid(this_z_indices,1:2);
x_grid = (1:n_grid)*(range(xy_locations(:,1))/n_grid) +min(xy_locations(:,1));
y_grid = (1:n_grid)*(range(xy_locations(:,2))/n_grid) +min(xy_locations(:,2));
% y_grid = (1:n_grid)*(80/n_grid) -40;
 cscale=   [min(true_shape_values) max(true_shape_values)];
[Xq,Yq] = meshgrid(x_grid,y_grid);
%             title(['z plane:' num2str(unique_z(i_z)) ])
%% Pick out the selected parameter sets
clear('param_list')
current_params=one_neuron.params(end);
%% Draw samples from the prior shape
[X_tmp,Y_tmp] = meshgrid(x_grid_sample,y_grid_sample);
x_vec=reshape(X_tmp,[],1);y_vec=reshape(Y_tmp,[],1);
new_locations= zeros(length(x_vec),3);
new_locations(:,1)=x_vec;new_locations(:,2)=y_vec;


%% Draw samples from each of the chosen iterations:
posterior_samples=cell([S 1]);
post_mean_shape=cell([1 1]);
clear('new_shape_params')
post_var_shape=cell([1 1]);
i=1;
[new_shape_params(i)]=get_new_shape_conditional(current_params,new_locations,prior_info);
for s=1:S
    posterior_sample = draw_samples_from_var_dist(current_params);
    post_new_sample=draw_samples_from_shape_conditional(new_shape_params(i),posterior_sample);
    posterior_samples{s,i}=[posterior_sample.shapes; post_new_sample.shapes];
    %=[diag(current_params.shapes.Sigma_tilde);diag(new_shape_params(i).Sigma_cond)];
end
post_mean_shape{i}= [posterior_samples{1,i}]/S;
for s=2:S
    post_mean_shape{i}= post_mean_shape{i}+[posterior_samples{s,i}]/S;
end
% calculate the posterior variance from the samples:
post_var_shape{i}=zeros(length(post_mean_shape{i}),1);

for i_loc = 1:length(post_mean_shape{i})
    tmp=zeros(S,1);
    for s=1:S
        tmp(s)=   posterior_samples{s,i}(i_loc);
    end
    post_var_shape{i}(i_loc)=var(tmp);
end

%% Posterior variances
Vq_fits = griddata(new_shape_params(end).locations(:,1),new_shape_params(end).locations(:,2),post_var_shape{end},Xq,Yq);
Vq_vec=[Vq_fits(:)];Xq_vec=[Xq(:)];Yq_vec=[Yq(:)];
optimal_locations=zeros(n_locations,3);

for i=1:n_locations 
    %
        [maxValue, max_index] = max(Vq_vec);
        optimal_locations(i,1:2)=[Xq_vec(max_index) Yq_vec(max_index)];
        neighbour_indices= find( abs(Xq_vec-Xq_vec(max_index))<gap &  abs(Yq_vec-Yq_vec(max_index))<gap   ) ;
        Vq_vec(neighbour_indices)=0;
    %
end
%% Contrasts
Vq = griddata(xy_locations(:,1),xy_locations(:,2),true_shape_values,Xq,Yq);
Vq_fits = griddata(new_shape_params(end).locations(:,1),new_shape_params(end).locations(:,2),post_mean_shape{end},Xq,Yq);
Vq_diff=abs(Vq_fits-Vq);
Vq_diff_vec=[Vq_diff(:)];
oracle_locations=zeros(n_locations,3);
for i=1:n_locations 
    %
        [maxValue, max_index] = max(Vq_diff_vec);
        oracle_locations(i,1:2)=[Xq_vec(max_index) Yq_vec(max_index)];
        neighbour_indices= find( abs(Xq_vec-Xq_vec(max_index))<gap &  abs(Yq_vec-Yq_vec(max_index))<gap   ) ;
        Vq_diff_vec(neighbour_indices)=0;
    %
end

