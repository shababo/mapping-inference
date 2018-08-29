function [GP_samples] = draw_1D_GP(locations,n_shapes,GP_params)
%% Draw 3D GP from GP_params
% GP_params should contain one model
% interpolation='linear';
% GP_params.mean_params .var_params .tau
interpolation='square';
epsilon=1e-3;
GP_samples=struct;
%% Find unique locations to reduce computing ost 
[locations_unique,~,ic] = unique(locations,'rows');
%%
switch interpolation
    case 'linear'
       
        
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

