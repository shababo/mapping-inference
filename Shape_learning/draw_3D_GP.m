function [GP_samples] = draw_3D_GP(locations,n_shapes, GP_params)
%% Draw 3D GP from GP_params
% GP_params should contain three models of the three dimensions
% axis_list={'x' 'y' 'z'};
axis_list=fieldnames(GP_params);
tmp_i=strcmp(axis_list,'type');
axis_list=axis_list(find(~tmp_i));
dims= size(locations,2);  % number of coordinates, 1 to 3
% interpolation='linear';
type=GP_params.type;
% 'xy_square' or 'square'
epsilon=1e-3;
GP_samples=struct;
% locations=simulation_params.mesh_grid;
% n_shapes=length(connected_cell_IDs);
% GP_params=prior_params.GP_params;

%% Find unique locations to reduce computing cost 
[locations_unique,~,ic] = unique(locations,'rows');
%%
[interpolated_shape]=interpolate_3D(locations_unique,GP_params,type);
% interpolated_shape$mean_3d, var_3d, GP_samples
%%
GP_samples = interpolated_shape.GP_samples;
GP_samples.full=struct;
GP_samples.full.locations=locations;
GP_samples.full.locations_unique=locations_unique;
GP_samples.full.mapping_unique=ic;

if ~strcmp(type,'square')
    % Draw samples from xy and z_separately
    
    GP_samples.full.Kcor_xy=GP_samples.x.Kcor.*GP_samples.y.Kcor;
    GP_samples.full.mean_xy=interpolated_shape.mean_xy;
    sigma_mat=(sqrt(interpolated_shape.var_xy)*ones(1,length( interpolated_shape.var_xy)));
    GP_samples.full.Kcov_xy= sigma_mat.*GP_samples.full.Kcor_xy.*sigma_mat';
    
    GP_samples.full.Kcor_z=GP_samples.z.Kcor;
    GP_samples.full.mean_z=interpolated_shape.mean_z;
    sigma_mat=(sqrt(interpolated_shape.var_z)*ones(1,length(interpolated_shape.var_z)));
    GP_samples.full.Kcov_z= sigma_mat.*GP_samples.full.Kcor_z.*sigma_mat';
    % Draw samples from xy and z:
    GP_samples.full.samples_xy = zeros(size(GP_samples.full.mean_xy, 1),n_shapes );
    GP_samples.full.samples_z = zeros(size(GP_samples.full.mean_z, 1),n_shapes );
    for i_shape = 1:n_shapes
        GP_samples.full.samples_xy(:,i_shape)=min(1,max(0,mvnrnd(GP_samples.full.mean_xy,GP_samples.full.Kcov_xy)))';
        GP_samples.full.samples_z(:,i_shape)=min(1,max(0,mvnrnd(GP_samples.full.mean_z,GP_samples.full.Kcov_z)))';
    end
    
    % Interpolate the 3D shapes from the 2D and 1D samples
    GP_samples.full.samples_unique = zeros(size(locations_unique,1),n_shapes );
      GP_samples.full.mean = zeros(size(locations_unique,1),1);GP_samples.full.var = zeros(size(locations_unique,1),1);
    for i_loc = 1:size(locations_unique,1)
        i_xy=interpolated_shape.index_xy(i_loc);
        i_z=interpolated_shape.index_z(i_loc);
        % Interpolate the xy-mean:
        sqloc=locations_unique(i_loc,:).^2;
        sqloc(sqloc==0)=epsilon;
        ssq = sum(sqloc');
        ssq_mat = ssq' *ones(1,3);
        exp_indices = sqloc./ssq_mat;
        exp_tmp = [1-exp_indices(3) exp_indices(3)];
 exp_mean=ones(1,2);
        for i_shape = 1:n_shapes
            sample_tmp=[GP_samples.full.samples_xy(i_xy,i_shape) GP_samples.full.samples_z(i_z,i_shape)];
            sample_tmp=(sample_tmp.^exp_mean);
            sample_3d=prod(sample_tmp')';
            GP_samples.full.samples_unique(i_loc,i_shape)=sample_3d;
            
            mean_tmp=[GP_samples.full.mean_xy(i_xy) GP_samples.full.mean_z(i_z)];
            mean_tmp=(mean_tmp.^exp_mean);
            mean_3d=prod(mean_tmp')';
            GP_samples.full.mean(i_loc)=mean_3d;
            
            var_tmp=[interpolated_shape.var_xy(i_xy) interpolated_shape.var_z(i_z)];
            var_tmp=(var_tmp.^exp_tmp);
            var_3d=prod(var_tmp')';
            GP_samples.full.var(i_loc)=var_3d;
        end
        
        %     GP_samples.full.loglklh(i_shape)= log(mvnpdf(GP_samples.full.samples_unique(:,i_shape),GP_samples.full.mean,GP_samples.full.Kcov));
    end
    
else
    for i_ax = 1:dims
        ax=axis_list{i_ax};
        if i_ax == 1
            GP_samples.full.Kcor=GP_samples.(ax).Kcor;
        else
            GP_samples.full.Kcor=GP_samples.full.Kcor.*GP_samples.(ax).Kcor;
        end
    end
    GP_samples.full.mean=interpolated_shape.mean_3d;
    sigma_mat=(sqrt(interpolated_shape.var_3d)*ones(1,length(interpolated_shape.var_3d)));
    GP_samples.full.Kcov= sigma_mat.*GP_samples.full.Kcor.*sigma_mat';
    % GP_samples.full.Kcov=(GP_samples.full.Kcov+GP_samples.full.Kcov')/2;
    GP_samples.full.samples_unique = zeros(size(locations_unique,1),n_shapes );
    GP_samples.full.loglklh=zeros(n_shapes,1 );
    for i_shape = 1:n_shapes
        GP_samples.full.samples_unique(:,i_shape)=min(1,max(0,mvnrnd(GP_samples.full.mean,GP_samples.full.Kcov)))';
        %     GP_samples.full.loglklh(i_shape)= log(mvnpdf(GP_samples.full.samples_unique(:,i_shape),GP_samples.full.mean,GP_samples.full.Kcov));
    end
    GP_samples.full.var=interpolated_shape.var_3d;
end

GP_samples.full.samples=GP_samples.full.samples_unique(ic,:);
%  Drawing shifts:
% GP_samples.full.shifts=zeros(dims,n_shapes);
% 
% for i_ax = 1:dims
%     ax = axis_list{i_ax};
%     switch GP_params.(ax).shift_params.type
%         case 'normal'
%             GP_samples.full.shifts(i_ax, :)=normrnd(GP_params.(ax).shift_params.mean,...
%                 GP_params.(ax).shift_params.var,[1 n_shapes]);
%     end
% end
