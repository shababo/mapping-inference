function [neurons] = assign_shapes_to_neurons(neurons,simulation_params,prior_info)
%%
n=length(neurons);
GP_params=prior_info.GP_params;
if simulation_params.batch.random_shape % false: use mean shape;
    tmp= draw_3D_GP(simulation_params.mesh_grid,n,GP_params);
    for i=1:n
        neurons(i).truth.shape= tmp.full.samples(:,i);
    end
else % use the mean shape
     tmp= draw_3D_GP(simulation_params.mesh_grid,1,GP_params);
    for i=1:n
        neurons(i).truth.shape= tmp.full.mean(tmp.full.mapping_unique);
    end
    
end
