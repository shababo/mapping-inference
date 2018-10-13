function [new_shape_sample] = draw_samples_from_shape_conditional(new_shape_params,posterior_sample)
%% Draw a sample from the NVM defined by parameters in the new_shape_params
%   .mean_new mean vector of the new shapes 
%   .mean_old mean vector of existing shapes
%   .Sigma_cross off-diagonal block (new old)
%   .Sigma_new covariance of the new shapes
%   .Sigma_old covariance of the old shapes
%   .Sigma_mean product of Simga_cross*inv(Sigma_old)
%   .Sigma_cond conditional covariance 
%   .locations % locaitons of the new shapes 
% We drop all locations that have been stimulated 


n_cell = length(posterior_sample);
clear('new_shape_sample')
new_shape_sample(n_cell)=struct;

for i_cell =1:n_cell
    this_param=new_shape_params(i_cell);
   old_shape_values = posterior_sample(i_cell).shapes; 
   new_shape_mean = this_param.mean_new - this_param.Sigma_mean*(old_shape_values- this_param.mean_old);
   
   new_shape_sample(i_cell).shapes=mvnrnd(new_shape_mean,this_param.Sigma_cond)';
   
end
