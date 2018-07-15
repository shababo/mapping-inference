function [mean_3d, var_3d, GP_samples] = interpolate_3D(locations,GP_params,type)
epsilon=1e-3;
axis_list=fieldnames(GP_params);
dims=length(axis_list);
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
            
        end
        sqloc(sqloc==0)=epsilon;
        ssq = sum(sqloc');
        ssq_mat = ssq' *ones(1,3);
        exp_indices = sqloc./ssq_mat;
        
        mean_tmp=(mean_val.^exp_indices);
        mean_3d=prod(mean_tmp')';
        
        var_tmp=(var_val.^exp_indices);
        var_3d=prod(var_tmp')';
        
end