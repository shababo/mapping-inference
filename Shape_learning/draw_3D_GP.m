function [GP_samples] = draw_3D_GP(locations,n_shapes, GP_params)
%% Draw 3D GP from GP_params
% GP_params should contain three models of the three dimensions
axis_list={'x' 'y' 'z'};
dims= size(locations,2);  % number of coordinates, 1 to 3


GP_samples=struct;
for i_ax = 1:dims
    ax = axis_list{i_ax};
    GP_samples.(ax)=struct;
    X=locations(:,i_ax);
    GP_samples.(ax).mean=quick_match(X,GP_params.(ax).mean_params);
    GP_samples.(ax).var=quick_match(X,GP_params.(ax).var_params);
    tau=GP_params.(ax).tau;
    
    GP_samples.(ax).Kcor=get_kernel_cor(X,X,tau);
    sigma_mat = sqrt(GP_samples.(ax).var)*ones(1,length(X));
    GP_samples.(ax).Kcov= sigma_mat.*GP_samples.(ax).Kcor.*sigma_mat';
    GP_samples.(ax).Kcov= (GP_samples.(ax).Kcov+GP_samples.(ax).Kcov')/2';
    
end

GP_samples.full=struct;
GP_samples.full.locations=locations;
for i_ax = 1:dims
   ax=axis_list{i_ax};
    if i_ax == 1
       
       GP_samples.full.mean=GP_samples.(ax).mean;
       GP_samples.full.var=GP_samples.(ax).var;
       GP_samples.full.Kcov=GP_samples.(ax).Kcov;
   else
       GP_samples.full.mean=GP_samples.full.mean.*GP_samples.(ax).mean;
       GP_samples.full.var=GP_samples.full.var.*GP_samples.(ax).var;
       GP_samples.full.Kcov=GP_samples.full.Kcov.*GP_samples.(ax).Kcov;
   end
end
GP_samples.full.samples = zeros(size(locations,1),n_shapes );
for i_shape = 1:n_shapes
    GP_samples.full.samples(:,i_shape)=mvnrnd(GP_samples.full.mean,GP_samples.full.Kcov)';
end

%% Drawing shifts:
GP_samples.full.shifts=zeros(dims,n_shapes);

for i_ax = 1:dims
    ax = axis_list{i_ax};
    switch GP_params.(ax).shift_params.type
        case 'normal'
            GP_samples.full.shifts(i_ax, :)=normrnd(GP_params.(ax).shift_params.mean,...
                GP_params.(ax).shift_params.var,[1 n_shapes]);
    end
end
