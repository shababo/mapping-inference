function [GP_samples] = draw_1D_GP(locations,n_shapes,GP_params)
%% Draw 3D GP from GP_params
% GP_params should contain three models of the three dimensions
% interpolation='linear';
interpolation='square';
epsilon=1e-3;
GP_samples=struct;
%% Find unique locations to reduce computing ost 
[locations_unique,~,ic] = unique(locations,'rows');
%%
switch interpolation
    case 'linear'
        mean_vec =[];
        loc_vec=zeros([0 3]);
        var_vec=[];
        range_increase=2; % double the min and max
        % we also want to increase the range of grid (fill with zeros) for
        % interpolation
        for i_ax = 1:dims
            ax = axis_list{i_ax};
            extended_grid =  GP_params.mean_params.grid;
            grid_gap = GP_params.mean_params.grid(2)-GP_params.mean_params.grid(1);
            right_ext= max(extended_grid):grid_gap:range_increase*max(extended_grid);
            left_ext = range_increase*min(extended_grid):grid_gap:min(extended_grid);
            extended_grid = [left_ext extended_grid right_ext];
            loc_tmp = zeros([length(extended_grid) 3]);
            
            loc_tmp(:,i_ax)= extended_grid;
            
            loc_vec=[loc_vec; loc_tmp];
            
            extended_mean=[zeros(length(left_ext),1);  GP_params.mean_params.values; zeros(length(right_ext),1)];
            mean_vec = [mean_vec; extended_mean];
            
            extended_var=[zeros(length(left_ext),1);  GP_params.var_params.values; zeros(length(right_ext),1)];
            var_vec = [var_vec;  extended_var];
            
            X=locations_unique(:,i_ax);
            
            tau=GP_params.tau;
            GP_samples.Kcor=get_kernel_cor(X,X,tau);
            GP_samples.Kcor= (GP_samples.Kcor+GP_samples.Kcor')/2;
            
        end
        %
        
        mean_3d = griddata(loc_vec(:,1),loc_vec(:,2),loc_vec(:,3), mean_vec,...
            locations_unique(:,1),locations_unique(:,2),locations_unique(:,3));
        
        var_3d = griddata(loc_vec(:,1),loc_vec(:,2),loc_vec(:,3), var_vec,...
           locations_unique(:,1),locations_unique(:,2),locations_unique(:,3));
        
    case 'square'
        % Interpolation as weighted averages
        
        mean_val=locations_unique;
        var_val=locations_unique;
        sqloc=locations_unique.^2;
        radius_to_center = sqrt(sqloc);
        
            
            X=locations_unique;
            tau=GP_params.tau;
            GP_samples.Kcor=get_kernel_cor(X,X,tau);
            GP_samples.Kcor= (GP_samples.Kcor+GP_samples.Kcor')/2;
            
            X_r=sign(X).*radius_to_center;
            mean_val=quick_match(X_r,GP_params.mean_params);
            var_val=quick_match(X_r,GP_params.var_params); 
            
        
end
%%
GP_samples.full=struct;
GP_samples.full.locations=locations;
GP_samples.full.locations_unique=locations_unique;
GP_samples.full.mapping_unique=ic;
GP_samples.full.Kcor=GP_samples.Kcor;
GP_samples.full.mean=mean_val;
sigma_mat=(sqrt(var_val)*ones(1,length(var_val)));
GP_samples.full.Kcov= sigma_mat.*GP_samples.full.Kcor.*sigma_mat';
GP_samples.full.samples_unique = zeros(size(locations_unique,1),n_shapes );
GP_samples.full.loglklh=zeros(n_shapes,1 );
for i_shape = 1:n_shapes
    GP_samples.full.samples_unique(:,i_shape)=max(0,mvnrnd(GP_samples.full.mean,GP_samples.full.Kcov))';
%     GP_samples.full.loglklh(i_shape)= log(mvnpdf(GP_samples.full.samples_unique(:,i_shape),GP_samples.full.mean,GP_samples.full.Kcov));
%     if isinf(GP_samples.full.loglklh(i_shape))
%     GP_samples.full.loglklh(i_shape)=10;
%     end
GP_samples.full.samples(:,i_shape)=GP_samples.full.samples_unique(ic,i_shape);

end

