function [interpolated_shape] = interpolate_3D(locations,GP_params,type)
%%
% locations=locations_unique;
%%
epsilon=1e-3;
axis_list=fieldnames(GP_params);

tmp_i=strcmp(axis_list,'type');
axis_list=axis_list(find(~tmp_i));
dims=min(3,length(axis_list));
switch type
    case 'linear'
        mean_vec =[];
        loc_vec=zeros([0 3]);
        var_vec=[];
        range_increase=2; % double the min and max
        % we also want to increase the range of grid (fill with zeros) for
        % interpolation
        for i_ax = 1:dims
            ax = axis_list{i_ax};
            extended_grid =  GP_params.(ax).mean_params.grid;
            grid_gap = GP_params.(ax).mean_params.grid(2)-GP_params.(ax).mean_params.grid(1);
            right_ext= max(extended_grid):grid_gap:range_increase*max(extended_grid);
            left_ext = range_increase*min(extended_grid):grid_gap:min(extended_grid);
            extended_grid = [left_ext extended_grid right_ext];
            loc_tmp = zeros([length(extended_grid) 3]);
            
            loc_tmp(:,i_ax)= extended_grid;
            
            loc_vec=[loc_vec; loc_tmp];
            
            extended_mean=[zeros(length(left_ext),1);  GP_params.(ax).mean_params.values; zeros(length(right_ext),1)];
            mean_vec = [mean_vec; extended_mean];
            
            extended_var=[zeros(length(left_ext),1);  GP_params.(ax).var_params.values; zeros(length(right_ext),1)];
            var_vec = [var_vec;  extended_var];
            
            X=locations(:,i_ax);
            
            tau=GP_params.(ax).tau;
            GP_samples.(ax).Kcor=get_kernel_cor(X,X,tau);
            GP_samples.(ax).Kcor= (GP_samples.(ax).Kcor+GP_samples.(ax).Kcor')/2;
            
        end
        %
        
        mean_3d = griddata(loc_vec(:,1),loc_vec(:,2),loc_vec(:,3), mean_vec,...
            locations_unique(:,1),locations_unique(:,2),locations_unique(:,3));
        
        var_3d = griddata(loc_vec(:,1),loc_vec(:,2),loc_vec(:,3), var_vec,...
            locations_unique(:,1),locations_unique(:,2),locations_unique(:,3));
        
    case 'square'
        % Interpolation as weighted averages
        
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
            mean_val(:,i_ax)=quick_match(X_r,GP_params.(ax).mean_params);
            var_val(:,i_ax)=quick_match(X_r,GP_params.(ax).var_params);
            %             mean_val(:,i_ax)=quick_match(X,GP_params.(ax).mean_params);
            %             var_val(:,i_ax)=quick_match(X,GP_params.(ax).var_params);
            
        end
        sqloc(sqloc==0)=epsilon;
        ssq = sum(sqloc');
        ssq_mat = ssq' *ones(1,3);
        exp_indices = sqloc./ssq_mat;
        
        mean_tmp=(mean_val.^exp_indices);
        mean_3d=prod(mean_tmp')';
        
        var_tmp=(var_val.^exp_indices);
        var_3d=prod(var_tmp')';
    case 'xy_square'
        % Interpolation as weighted averages
        xy_locations = locations(:,1:2);
        z_locations = locations(:,3);
        [xy_unique,~,i_xy] = unique(xy_locations,'rows');
        [z_unique,~,i_z] = unique(z_locations,'rows');
        sqloc=xy_unique.^2;
        radius_to_center = sqrt(sum(sqloc'))';
        mean_val=struct;var_val=struct;
        for i_ax = 1:dims
            ax = axis_list{i_ax};
            if i_ax ==3
                X=z_unique;
            else
                X=xy_unique(:,i_ax);
                %                 X=sign(X).*radius_to_center;
            end
            tau=GP_params.(ax).tau;
            GP_samples.(ax).Kcor=get_kernel_cor(X,X,tau);
            GP_samples.(ax).Kcor= (GP_samples.(ax).Kcor+GP_samples.(ax).Kcor')/2;
            
            mean_val.(ax)=quick_match(X,GP_params.(ax).mean_params);
            var_val.(ax)=quick_match(X,GP_params.(ax).var_params);
        end
        
        % Interpolate the xy-mean:
        sqloc=xy_unique.^2;
        sqloc(sqloc==0)=epsilon;
        ssq = sum(sqloc');
        ssq_mat = ssq' *ones(1,2);
        exp_indices = sqloc./ssq_mat;
        exp_mean=ones(1,2);
        mean_tmp=[mean_val.x mean_val.y];
        mean_tmp=(mean_tmp.^exp_mean);
        xy_mean=prod(mean_tmp')';
        
        var_tmp=[var_val.x var_val.y];
        var_tmp=(var_tmp.^exp_indices);
        xy_var=prod(var_tmp')';
        
    case 'xy'
        %
    xy_locations = locations(:,1:2);
        z_locations = locations(:,3);
        [xy_unique,~,i_xy] = unique(xy_locations,'rows');
        [z_unique,~,i_z] = unique(z_locations,'rows');
        
        sqloc=xy_unique.^2;
        radius_to_center = sqrt(sum(sqloc'))';
        mean_val=struct;var_val=struct;
        for i_ax = 1:3
            ax = axis_list{i_ax};
            if i_ax ==3
                X=z_unique;
                tau=GP_params.(ax).tau;
            else
                X=xy_unique(:,i_ax);
                tau=GP_params.xy.tau(i_ax);
            end
            GP_samples.(ax).Kcor=get_kernel_cor(X,X,tau);
            GP_samples.(ax).Kcor= (GP_samples.(ax).Kcor+GP_samples.(ax).Kcor')/2;
        end
        mean_val.z=quick_match(z_unique,GP_params.z.mean_params);
        var_val.z=quick_match(z_unique,GP_params.z.var_params);
        
        xy_mean=quick_match(xy_unique,GP_params.xy.mean_params);
        xy_var=quick_match(xy_unique,GP_params.xy.var_params);
        
        
end
%% Save output
interpolated_shape = struct;
interpolated_shape.type = type;
if ~strcmp(type,'square')
    interpolated_shape.locations = locations;
    
    interpolated_shape.mean_xy = xy_mean;
    interpolated_shape.var_xy = xy_var;
    interpolated_shape.index_xy = i_xy;
    
    interpolated_shape.mean_z =  mean_val.z;
    interpolated_shape.var_z =  var_val.z;
    interpolated_shape.index_z = i_z;
    
    interpolated_shape.GP_samples = GP_samples;
else
    interpolated_shape.locations = locations;
    interpolated_shape.mean_3d = mean_3d;
    interpolated_shape.var_3d = var_3d;
    interpolated_shape.GP_samples = GP_samples;
end
