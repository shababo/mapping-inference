function [new_shape_params]=get_new_shape_conditional(neurons,new_trials,prior_info)
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
n_cell=length(neurons);
clear('new_shape_params');
new_shape_params(n_cell)=struct;
for i=1:n_cell
new_shape_params(i).locations=zeros(0,3);
end
n_new=zeros(n_cell,1);
boundary_params = prior_info.prior_parameters.boundary_params;
for i_trial = 1:length(new_trials)
    stim_locs = new_trials(i_trial).locations;
    for i_loc = size(stim_locs,1)
        if ~isnan(stim_locs(i_loc,1))
            for i_cell = 1:n_cell
                this_loc=neurons(i_cell).location;
                rel_pos=stim_locs(i_loc,:)-this_loc;
                if check_in_boundary(rel_pos,boundary_params)
                    [C,ia,ib] = intersect(neurons(i_cell).params(end).shapes.locations,rel_pos,'rows');
                    if isempty(C) % this is a new location:
                        [C,ia,ib] = intersect(new_shape_params(i_cell).locations,rel_pos,'rows');
                    if isempty(C) % this is the first time we see it in new_trials:
                        n_new(i_cell)=n_new(i_cell)+1;
                        new_shape_params(i_cell).locations(n_new(i_cell),:)=rel_pos;
                    end     % do nothing
                    end
                end
                
            end
        end
    end
end

%% Calculate the parameters in the conditional distribution (correlation matrix for the shapes) 

for i_cell = 1:n_cell
    if ~isempty(  new_shape_params(i_cell).locations)
        locs_old=neurons(i_cell).params(end).shapes.locations;
        locs_new=  new_shape_params(i_cell).locations;
        locs=[locs_new; locs_old];
        new_old_flags = [ones(size(locs_new,1),1); zeros(size(locs_old,1),1)];
        [mean_3d, var_3d]=interpolate_3D(locs,GP_params,type);
        new_shape_params(i_cell).mean_new=mean_3d(new_old_flags==1);
        new_shape_params(i_cell).mean_old=mean_3d(new_old_flags==0);
        
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
        Full_Kcov= sigma_mat.*Full_Kcor.*sigma_mat';
        
        new_shape_params(i_cell).Sigma_cross= Full_Kcov(new_old_flags==1,new_old_flags==0);
        new_shape_params(i_cell).Sigma_new= Full_Kcov(new_old_flags==1,new_old_flags==1);
         new_shape_params(i_cell).Sigma_old= Full_Kcov(new_old_flags==0,new_old_flags==0);
        inv_tmp=inv(new_shape_params(i_cell).Sigma_old);
         new_shape_params(i_cell).Sigma_mean= new_shape_params(i_cell).Sigma_cross*inv_tmp;
        
         new_shape_params(i_cell).Sigma_cond= new_shape_params(i_cell).Sigma_new-...
             new_shape_params(i_cell).Sigma_cross*inv_tmp*new_shape_params(i_cell).Sigma_cross';
    end
end