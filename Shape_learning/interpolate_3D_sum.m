function [val_3d]=interpolate_3D_sum(locations, xy_val,i_xy,z_val,i_z,epsilon)
 % Interpolate the 3D GP values: 
        tmp=locations.^2;
        sqloc=[tmp(1)+tmp(2) tmp(3)];
        sqloc(sqloc==0)=epsilon;
        
        ssq = sum(sqloc');
        ssq_mat = ssq' *ones(1,2);
        exp_indices = sqloc./ssq_mat;
        
        xy_val_full = xy_val(i_xy);z_val_full =z_val(i_z); 
        val_tmp=([xy_val_full z_val_full].*exp_indices);
        val_3d=sum(val_tmp')';
        