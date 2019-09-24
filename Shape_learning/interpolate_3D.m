function [interpolated_shape] = interpolate_3D(locations,GP_params,type)
%%
% locations=locations_unique;
%%
epsilon=1e-3;
axis_list=fieldnames(GP_params);

tmp_i=strcmp(axis_list,'type');
dims=3;
switch type
    case 'linear'
        
    case 'x-y-z' % used to be square
        % Interpolate 3d mean & var as weighted product of three GPs
        axis_list={'x','y','z'};
        mean_val=locations;
        var_val=locations;
        sqloc=locations.^2;
        radius_to_center = sqrt(sum(sqloc'))';
        for i_ax = 1:dims
            ax = axis_list{i_ax};
            X=locations(:,i_ax);
            tau=GP_params.(ax).tau;
            GP_samples.(ax).Kcor=get_kernel_cor(X,X,tau);
            GP_samples.(ax).Kcor= (GP_samples.(ax).Kcor+GP_samples.(ax).Kcor')/2;
            X_r=sign(X).*radius_to_center;
            mean_val(:,i_ax)=GP_params.(ax).mean_params.excite(quick_match(X_r,GP_params.(ax).mean_params));
            var_val(:,i_ax)=quick_match(X_r,GP_params.(ax).var_params);
        end
        sqloc(sqloc==0)=epsilon;
        ssq = sum(sqloc');
        ssq_mat = ssq' *ones(1,3);
        exp_indices = sqloc./ssq_mat;
        mean_tmp=(mean_val.^exp_indices);
        mean_3d=prod(mean_tmp')';
        var_tmp=(var_val.^exp_indices);
        var_3d=prod(var_tmp')';
    case '(x-y)-z' % used to be xy-square
        % Interpolate 3d mean & var as weighted product of three GPs
        % Output the interpolation of xy and z
        axis_list={'x','y','z'};
        mean_val=struct;var_val=struct;
        xy_locations = locations(:,1:2);
        z_locations = locations(:,3);
        [xy_unique,~,i_xy] = unique(xy_locations,'rows');
        [z_unique,~,i_z] = unique(z_locations,'rows');
        
        for i_ax = 1:dims
            ax = axis_list{i_ax};
            if i_ax == 3
            X=xy_unique(:,i_ax);
            else
            X=z_unique;    
            end
            tau=GP_params.(ax).tau;
            mean_val.(ax)=GP_params.(ax).mean_params.excite(quick_match(X,GP_params.(ax).mean_params));
            var_val.(ax)=quick_match(X,GP_params.(ax).var_params);
            if i_ax == 1
                GP_samples.xy.Kcor=get_kernel_cor(X,X,tau);
                
            elseif i_ax == 2
                GP_samples.xy.Kcor=GP_samples.xy.Kcor.*get_kernel_cor(X,X,tau);
                GP_samples.xy.Kcor= (GP_samples.xy.Kcor+GP_samples.xy.Kcor')/2;
            end
        end
        GP_samples.z.Kcor=get_kernel_cor(X,X,tau);
        GP_samples.z.Kcor= (GP_samples.z.Kcor+GP_samples.z.Kcor')/2;
        
   % Interpolate the xy-mean:
        sqloc=locations.^2;
        sqloc(sqloc==0)=epsilon;
        ssq = sum(sqloc(1:2)'); ssq_mat = ssq' *ones(1,2);
        exp_indices = sqloc(1:2)./ssq_mat;
        
        mean_tmp=[mean_val.x mean_val.y];
        mean_tmp=(mean_tmp.^exp_indices);
        xy_mean=prod(mean_tmp')';
        var_tmp=[var_val.x var_val.y];
        var_tmp=(var_tmp.^exp_indices);
        xy_var=prod(var_tmp')';
        
        
        % Interpolate the 3D GP values:
        ssq = sum(sqloc'); ssq_mat = ssq' *ones(1,2);
        exp_indices = [sqloc(1)+sqloc(2) sqloc(3)]./ssq_mat;
        

        
        xy_mean_full = xy_mean(i_xy);z_mean_full = z_mean(i_z); 
        mean_tmp=([xy_mean_full z_mean_full].^exp_indices);
        mean_3d=prod(mean_tmp')';
        
        
        xy_var_full = xy_var(i_xy);z_var_full = z_var(i_z); 
        var_tmp=([xy_var_full z_var_full].^exp_indices);
        var_3d=prod(var_tmp')';
     
    case 'xy-z' % used to be xy
        % Obtain the values from the xy GP and z GP
        % Interpolate the 3D value
        
        xy_locations = locations(:,1:2);
        z_locations = locations(:,3);
        [xy_unique,~,i_xy] = unique(xy_locations,'rows');
        [z_unique,~,i_z] = unique(z_locations,'rows');
        
        z_mean=quick_match(z_unique,GP_params.z.mean_params);
        z_var=quick_match(z_unique,GP_params.z.var_params);
        xy_mean=quick_match(xy_unique,GP_params.xy.mean_params);
        xy_var=quick_match(xy_unique,GP_params.xy.var_params);
       
        x=xy_unique(:,1);y=xy_unique(:,2);z=z_unique;
        GP_samples.xy.Kcor=get_kernel_cor(x,x,GP_params.xy.tau(1)).*get_kernel_cor(y,y,GP_params.xy.tau(2));
        GP_samples.xy.Kcor= (GP_samples.xy.Kcor+GP_samples.xy.Kcor')/2;    
        GP_samples.z.Kcor=get_kernel_cor(z,z,GP_params.z.tau);
        GP_samples.z.Kcor= (GP_samples.z.Kcor+GP_samples.z.Kcor')/2;

       [mean_3d]=GP_params.interpolation_func(locations, xy_mean,i_xy,z_mean,i_z,epsilon);
       
       [var_3d]=GP_params.interpolation_func(locations, xy_var,i_xy,z_var,i_z,epsilon);
end
%% Save output
interpolated_shape = struct;
interpolated_shape.type = type;
if strcmp(type,'x-y-z')
    interpolated_shape.locations = locations;
    interpolated_shape.mean_3d = mean_3d;
    interpolated_shape.var_3d = var_3d;
    interpolated_shape.GP_samples = GP_samples;
else
    
    interpolated_shape.locations = locations;
    
    interpolated_shape.mean_xy = xy_mean;
    interpolated_shape.var_xy = xy_var;
    interpolated_shape.index_xy = i_xy;
    interpolated_shape.mean_z =  z_mean;
    interpolated_shape.var_z =  z_var;
    interpolated_shape.index_z = i_z;
    interpolated_shape.mean_3d = mean_3d;
    interpolated_shape.var_3d =  var_3d;
    interpolated_shape.GP_samples = GP_samples;
end
