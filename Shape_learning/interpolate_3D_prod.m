function interpolate_3D_prod(locations, xy_mean,z_mean,epsilon)
 % Interpolate the 3D GP values: 
        tmp=locations.^2;
        sqloc=[tmp(1)+tmp(2) tmp(3)];
        sqloc(sqloc==0)=epsilon;
        
        ssq = sum(sqloc');
        ssq_mat = ssq' *ones(1,2);
        exp_indices = sqloc./ssq_mat;
        
        xy_mean_full = xy_mean(i_xy);z_mean_full =z_mean(i_z); 
        mean_tmp=([xy_mean_full z_mean_full].^exp_indices);
        mean_3d=prod(mean_tmp')';
        
        xy_var_full = xy_var(i_xy);z_var_full = z_var(i_z); 
        var_tmp=([xy_var_full z_var_full].^exp_indices);
        var_3d=prod(var_tmp')';