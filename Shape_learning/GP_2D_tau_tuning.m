function [GP_samples] = GP_2D_tau_tuning(taus, locations,xy_pilot)

mean_val=quick_match2D(locations,xy_pilot.mean_params);
var_val=quick_match2D(locations,xy_pilot.var_params);

sigma_mat = sqrt(var_val)*ones(1,length(mean_val));
K_mat=cell([2 1]);
for i=1:2
    K_mat{i}=get_kernel_cor(locations(:,i),locations(:,i),taus(i));
end
Kcov=sigma_mat.*K_mat{1}.*K_mat{2}.*sigma_mat';

this_sample=max(0,mvnrnd(mean_val,Kcov))';
GP_samples=struct;
GP_samples.sample=this_sample;
GP_samples.locations=locations;
GP_samples.taus=taus;