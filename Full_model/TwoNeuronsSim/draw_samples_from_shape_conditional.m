function [new_shape_sample] = draw_samples_from_shape_conditional(new_shape_params,posterior_sample)
%% Draw a sample from the NVM defined by parameters in the new_shape_params
%   .mean_new mean vector of the new shapes 
%   .mean_old mean vector of existing shapes
%   .Sigma_cross off-diagonal block (new old)
%   .Sigma_new covariance of the new shapes
%   .Sigma_old covariance of the old shapes
%   .Sigma_mean product of Simga_cross*inv(Sigma_old)
%   .Sigma_cond conditional covariance 
%   .locations % locations of the new shapes 

% Note that the Gaussian r.v.s are unbounded
% We apply a recifier to force the shape value to live between 0 and 1.

n_cell = length(posterior_sample);
clear('new_shape_sample')
new_shape_sample(n_cell)=struct;


%if strcmp(prior_info.GP_params.type, 'xy_square')
    n_axis_model =2;
    axis_names = {'xy', 'z'};
    % Calculate the joint distribution (correlation matrix for the shapes)
    if strcmp(new_shape_params(1).(axis_names{1}).dist,'mvn')
        for i_cell = 1:n_cell
            this_param=new_shape_params(i_cell);
            for i_axis = 1:n_axis_model
                if  ~isempty(this_param.(axis_names{i_axis}).mean_new)
                old_shape_values = posterior_sample(i_cell).(axis_names{i_axis});
                new_shape_mean = this_param.(axis_names{i_axis}).mean_new + this_param.(axis_names{i_axis}).Sigma_mean*(old_shape_values- this_param.(axis_names{i_axis}).mean_old); %
                tmp=mvnrnd(new_shape_mean,this_param.(axis_names{i_axis}).Sigma_cond)';
                tmp(tmp<0)=0;tmp(tmp>1)=1;
                new_shape_sample(i_cell).(axis_names{i_axis})= tmp;
                else
                    new_shape_sample(i_cell).(axis_names{i_axis})=[];
                end
            end
        end
    end
%end
