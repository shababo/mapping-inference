function [new_shape_params]=get_new_shape_conditional(current_params,new_shape_params,prior_info)
%% Output: new_shape_params
%   .mean_new mean vector of the new shapes
%   .mean_old mean vector of existing shapes
%   .Sigma_cross off-diagonal block (new old)
%   .Sigma_new covariance of the new shapes
%   .Sigma_old covariance of the old shapes
%   .Sigma_mean product of Simga_cross*inv(Sigma_old)
%   .Sigma_cond conditional covariance
%   .locations % locaitons of the new shapes
%% Calculate the parameters in the conditional distribution (correlation matrix for the shapes)
GP_params=prior_info.GP_params;
boundary_params = GP_params.GP_boundary;
type=GP_params.type;
n_cell=length(current_params);
if strcmp(type, 'xy_square')
    n_axis_model =2;
    axis_names = {'xy', 'z'};
    
    % Calculate the joint distribution (correlation matrix for the shapes)
    
    for i_cell = 1:n_cell
        for i_axis = 1:n_axis_model
            if strcmp(current_params(1).(axis_names{i_axis}).dist,'mvn')
                locs_old=current_params(i_cell).(axis_names{i_axis}).locations;
                locs_new=  new_shape_params(i_cell).(axis_names{i_axis}).locations;
                locs=[locs_old;locs_new];
                
                new_shape_params(i_cell).(axis_names{i_axis}).all_locations=locs;
                new_shape_params(i_cell).(axis_names{i_axis}).old_locations=locs_old;
                new_old_flags = [zeros(size(locs_old,1),1); ones(size(locs_new,1),1)];
                
                % Fill the other locs with zeros before calling the
                % interpolate function: 
                
                if i_axis == 1
                    locs_3D=[locs ones(size(locs,1) ,1)];
                else
                    locs_3D=[ones(size(locs,1) , 2) locs ];
                end
                [interpolated_shape]=interpolate_3D(locs_3D,GP_params,type);
                
  
                if i_axis == 1
                    mean_3d=interpolated_shape.mean_xy;
                    var_3d=interpolated_shape.var_xy;
                else
                    mean_3d=interpolated_shape.mean_z;
                    var_3d=interpolated_shape.var_z;
                end
                 new_shape_params(i_cell).(axis_names{i_axis}).mean_new=mean_3d(new_old_flags==1);
                new_shape_params(i_cell).(axis_names{i_axis}).mean_old=mean_3d(new_old_flags==0);
              
                % Obtain the covariance matrix
                X=locs;
                tau=GP_params.(axis_names{i_axis}).tau;
                if i_axis == 1
                    for i_tmp =1:2
                        tmp_Kcor=get_kernel_cor(X(:,i_tmp),X(:,i_tmp),tau(i_tmp));
                        tmp_Kcor=(tmp_Kcor+tmp_Kcor')/2;
                        if i_tmp == 1
                            Full_Kcor=tmp_Kcor;
                        else
                            Full_Kcor=Full_Kcor.*tmp_Kcor;
                        end
                    end
                else
                    tmp_Kcor=get_kernel_cor(X,X,tau);
                    tmp_Kcor=(tmp_Kcor+tmp_Kcor')/2;
                    Full_Kcor=tmp_Kcor;
                end
                
                
                sigma_mat=sqrt(var_3d)*ones(1,length(var_3d));
                Full_Kcov= sigma_mat.*Full_Kcor.*sigma_mat';
                
                new_shape_params(i_cell).(axis_names{i_axis}).Sigma_cross= Full_Kcov(new_old_flags==1,new_old_flags==0);
                new_shape_params(i_cell).(axis_names{i_axis}).Sigma_new= Full_Kcov(new_old_flags==1,new_old_flags==1);
                new_shape_params(i_cell).(axis_names{i_axis}).Sigma_old= Full_Kcov(new_old_flags==0,new_old_flags==0);
                inv_tmp=inv( new_shape_params(i_cell).(axis_names{i_axis}).Sigma_old);
                new_shape_params(i_cell).(axis_names{i_axis}).Sigma_mean=  new_shape_params(i_cell).(axis_names{i_axis}).Sigma_cross*inv_tmp;
                
                new_shape_params(i_cell).(axis_names{i_axis}).Sigma_cond=  new_shape_params(i_cell).(axis_names{i_axis}).Sigma_new-...
                    new_shape_params(i_cell).(axis_names{i_axis}).Sigma_cross*inv_tmp* new_shape_params(i_cell).(axis_names{i_axis}).Sigma_cross';
                new_shape_params(i_cell).(axis_names{i_axis}).Sigma_cond=( new_shape_params(i_cell).(axis_names{i_axis}).Sigma_cond+ new_shape_params(i_cell).(axis_names{i_axis}).Sigma_cond')/2;
                
            end
            
        end
    end
end