function [GP_samples] = draw_3D_GP(locations,n_shapes, GP_params)
%% Draw 3D GP from GP_params
% GP_params should contain three models of the three dimensions
% axis_list={'x' 'y' 'z'};
axis_list=fieldnames(GP_params);
dims= size(locations,2);  % number of coordinates, 1 to 3
% interpolation='linear';
type='square';
epsilon=1e-3;

GP_samples=struct;

%% Find unique locations to reduce computing ost 

[locations_unique,~,ic] = unique(locations,'rows');

%%
[mean_3d, var_3d, GP_samples]=interpolate_3D(locations,GP_params,type);

%%
GP_samples.full=struct;
GP_samples.full.locations=locations;
GP_samples.full.locations_unique=locations_unique;
GP_samples.full.mapping_unique=ic;

for i_ax = 1:dims
    ax=axis_list{i_ax};
    if i_ax == 1
        
        GP_samples.full.Kcor=GP_samples.(ax).Kcor;
    else
        GP_samples.full.Kcor=GP_samples.full.Kcor.*GP_samples.(ax).Kcor;
    end
end
GP_samples.full.mean=mean_3d;
sigma_mat=(sqrt(var_3d)*ones(1,length(var_3d)));
GP_samples.full.Kcov= sigma_mat.*GP_samples.full.Kcor.*sigma_mat';
% GP_samples.full.Kcov=(GP_samples.full.Kcov+GP_samples.full.Kcov')/2;
GP_samples.full.samples_unique = zeros(size(locations_unique,1),n_shapes );
GP_samples.full.loglklh=zeros(n_shapes,1 );
for i_shape = 1:n_shapes
    GP_samples.full.samples_unique(:,i_shape)=max(0,mvnrnd(GP_samples.full.mean,GP_samples.full.Kcov))';
%     GP_samples.full.loglklh(i_shape)= log(mvnpdf(GP_samples.full.samples_unique(:,i_shape),GP_samples.full.mean,GP_samples.full.Kcov));
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
