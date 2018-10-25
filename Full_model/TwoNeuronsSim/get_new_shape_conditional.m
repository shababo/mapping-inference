function [new_shape_params]=get_new_shape_conditional(parameter_current,new_locations,prior_info)
%% Output: new_shape_params
%   .mean_new mean vector of the new shapes
%   .mean_old mean vector of existing shapes
%   .Sigma_cross off-diagonal block (new old)
%   .Sigma_new covariance of the new shapes
%   .Sigma_old covariance of the old shapes
%   .Sigma_mean product of Simga_cross*inv(Sigma_old)
%   .Sigma_cond conditional covariance
%   .locations % locaitons of the new shapes
% We drop all locations that have been stimulated
type='square';
GP_params=prior_info.prior_parameters.GP_params;
%% Detect new locations:
clear('new_shape_params');
min_gap = 2; % 2 um apart!
sigma_epsilon =  0.1; % a small value that mimics the noise 
new_shape_params=struct;
new_shape_params.locations=zeros(0,3);
n_new=0;
boundary_params = prior_info.prior_parameters.boundary_params;

this_loc=[0 0 0];
locs_old=parameter_current.shapes.locations;
for i_loc = 1:size(new_locations,1)
    if ~isnan(new_locations(i_loc,1))
        rel_pos=new_locations(i_loc,:)-this_loc;
        if check_in_boundary(rel_pos,boundary_params)
            [C,ia,ib] = intersect(locs_old,rel_pos,'rows');
            min_dist= min(sum([(rel_pos-locs_old).^2]'));
            if isempty(C) & min_dist>min_gap % this is a new location:
                [C,ia,ib] = intersect(locs_old,rel_pos,'rows');
                if isempty(C) % this is the first time we see it in new_trials:
                    n_new=n_new+1;
                    new_shape_params.new_locations(n_new,:)=rel_pos;
                end     % do nothing
            end
        end
        
    end
end

%% Calculate the parameters in the conditional distribution (correlation matrix for the shapes)

if ~isempty(  new_shape_params.new_locations)
    locs_old=parameter_current.shapes.locations;
    locs_new=  new_shape_params.new_locations;
    locs=[locs_old;locs_new];
    new_shape_params.locations=locs;
    
    new_shape_params.old_locations=locs_old;
    new_old_flags = [zeros(size(locs_old,1),1); ones(size(locs_new,1),1)];
    [mean_3d, var_3d]=interpolate_3D(locs,GP_params,type);
    new_shape_params.mean_new=mean_3d(new_old_flags==1);
    new_shape_params.mean_old=mean_3d(new_old_flags==0);
    
    % Obtain the covariance matrix
    axis_list= fieldnames(GP_params);
    for i_ax = 1:length(axis_list)
        X=locs(:,i_ax);
        ax=axis_list{i_ax};
        tau=GP_params.(ax).tau;
        tmp_Kcor=get_kernel_cor(X,X,tau);
        tmp_Kcor=(tmp_Kcor+tmp_Kcor')/2;
        if i_ax == 1
            Full_Kcor=tmp_Kcor;
        else
            Full_Kcor=Full_Kcor.*tmp_Kcor;
        end
    end
    sigma_mat=sqrt(var_3d)*ones(1,length(var_3d));
    Full_Kcov= sigma_mat.*(Full_Kcor+sigma_epsilon*diag(ones(length(var_3d),1)) ).*sigma_mat';
    
    new_shape_params.Sigma_cross= Full_Kcov(new_old_flags==1,new_old_flags==0);
    new_shape_params.Sigma_new= Full_Kcov(new_old_flags==1,new_old_flags==1);
    new_shape_params.Sigma_old= Full_Kcov(new_old_flags==0,new_old_flags==0);
    inv_tmp=inv(new_shape_params.Sigma_old);
    new_shape_params.Sigma_mean= new_shape_params.Sigma_cross*inv_tmp;
    
    new_shape_params.Sigma_cond= new_shape_params.Sigma_new-...
        new_shape_params.Sigma_cross*inv_tmp*new_shape_params.Sigma_cross';
    new_shape_params.Sigma_cond=(new_shape_params.Sigma_cond+new_shape_params.Sigma_cond')/2;
end