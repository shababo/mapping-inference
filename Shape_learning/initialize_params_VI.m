function [variational_params, prior_params,trials]=initialize_params_VI(neurons,trials,inference_params,prior_info)
%% This function needs to be cleaned up
n_cell=length(neurons);
clear('variational_params')
for i_cell = 1:n_cell
    variational_params(i_cell)=neurons(i_cell).params(end);
end
zero_bound = 1e-10;

switch inference_params.shape_type
    case 'fact'
        
    case 'non-fact'
        % Non-factorize shape from the same prior:
        GP_params=prior_info.GP_params;
        type=GP_params.type;
        min_dist = GP_params.min_dist; % minimun distance to distinguish two stim spot
        boundary_params = prior_info.GP_params.GP_boundary;
        
        axis_names = {'xy', 'z'};
        n_axis_model=length(axis_names);
        for i_trial = 1:length(trials)
            stim_locs = trials(i_trial).locations;
            for i_loc = size(stim_locs,1)
                if ~isnan(stim_locs(i_loc,1))
                    trials(i_trial).cell_and_pos{i_loc}= zeros(0,2); %i_cell, i_pos
                    for i_cell = 1:n_cell
                        this_loc=neurons(i_cell).location;
                        rel_pos=stim_locs(i_loc,:)-this_loc;
                        if check_in_boundary(rel_pos,boundary_params)
                            % Check for xy shape and z shape:
                            existing_loc_dim=size(variational_params(i_cell).shapes.locations);
                            sq_dist=sum((variational_params(i_cell).shapes.locations- ones(existing_loc_dim(1),1)*rel_pos).^2,2);
                            if (existing_loc_dim(1) == 0) | min(sq_dist.^(1/2))> min_dist
                                tmp = [];
                            else
                                [~,tmp]=min(sq_dist.^(1/2));
                            end
                            %[C,tmp,~] = intersect(variational_params(i_cell).(axis_names{i_axis}).locations,this_rel_pos,'rows');
                            if isempty(tmp) % this is a new location:
                                variational_params(i_cell).shapes.locations=...
                                    [variational_params(i_cell).shapes.locations; rel_pos];
                                
                                [interpolated_shape]=interpolate_3D(rel_pos,GP_params,type);
                                
                                mean_3d=interpolated_shape.mean_3d;
                                var_3d=interpolated_shape.var_3d;
                                
                                % Adding a small variance to allow the shape to
                                % learn from the data set when the prior shape
                                % variance is too small.
                                var_3d_vi=var_3d;
%                                  if sqrt(sum(rel_pos.^2))>min_dist % if this loc is not near the nucleus
                                    if (var_3d <  prior_info.GP_params.GP_minimal_variance)
                                        var_3d =   prior_info.GP_params.GP_minimal_variance;
                                        var_3d_vi=var_3d;
                                    end
                                    if isfield(prior_info.GP_params, 'GP_added_variance')
                                        var_3d_vi = var_3d+prior_info.GP_params.GP_added_variance;
                                    end
%                                  end
                                
                                lower_bound =mean_3d-   inference_params.GP_shape_bound_factor*sqrt(var_3d_vi);
                                upper_bound =mean_3d+   inference_params.GP_shape_bound_factor*sqrt(var_3d_vi);
%                                  lower_bound =-10;upper_bound =0; % for debugging
                                variational_params(i_cell).shapes.bounds.low = [variational_params(i_cell).shapes.bounds.low;lower_bound];
                                variational_params(i_cell).shapes.bounds.up = [variational_params(i_cell).shapes.bounds.up;upper_bound];
                                % logit transform:
                                switch  variational_params(i_cell).shapes.dist
                                    case 'mvn'
                                        variational_params(i_cell).shapes.prior_sigma=[variational_params(i_cell).shapes.prior_sigma; ...
                                            var_3d];
                                        variational_params(i_cell).shapes.vi_sigma=[variational_params(i_cell).shapes.vi_sigma; ...
                                            var_3d_vi];
                                        variational_params(i_cell).shapes.mean=[variational_params(i_cell).shapes.mean; ...
                                            mean_3d];
                                        variational_params(i_cell).shapes.log_sigma=log(variational_params(i_cell).shapes.vi_sigma);
                                    case 'mvn-logit'
                                        variational_params(i_cell).shapes.prior_sigma=[variational_params(i_cell).shapes.prior_sigma; ...
                                            var_3d];
                                        variational_params(i_cell).shapes.vi_sigma=[variational_params(i_cell).shapes.vi_sigma; ...
                                            var_3d_vi];
                                        variational_params(i_cell).shapes.mean=[variational_params(i_cell).shapes.mean; ...
                                            log( (mean_3d-lower_bound)/(upper_bound-mean_3d))];
                                        variational_params(i_cell).shapes.log_sigma=log(variational_params(i_cell).shapes.vi_sigma);
                                end
                                ia=length(variational_params(i_cell).shapes.mean);
                            else
                                ia=tmp;
                            end
                        trials(i_trial).cell_and_pos{i_loc}=[trials(i_trial).cell_and_pos{i_loc};...
                            i_cell ia];
                        
                        end
                    end
                end
            end
        end
        % Calculate the joint distribution (correlation matrix for the shapes)
        if strcmp(variational_params(1).shapes.dist,'mvn')
            for i_cell = 1:n_cell
                
               X=variational_params(i_cell).shapes.locations;
                for i_axis = 1:n_axis_model
                    tau=GP_params.(axis_names{i_axis}).tau;
                    if i_axis == 1
                        for i_tmp =1:2
                            tmp_Kcor=get_kernel_cor(X(:,i_tmp),X(:,i_tmp),tau(i_tmp));
                            tmp_Kcor=(tmp_Kcor+tmp_Kcor')/2;
                            if i_tmp == 1
                                Full_Kcor=tmp_Kcor;
                            else
                                Full_Kcor=Full_Kcor.*tmp_Kcor;
                            end
                        end
                    else
                        tmp_Kcor=get_kernel_cor(X(:,3),X(:,3),tau);
                        tmp_Kcor=(tmp_Kcor+tmp_Kcor')/2;
                        Full_Kcor=Full_Kcor.*tmp_Kcor;
                    end
                end
                Full_Kcor=GP_params.nugget*diag(ones([size(X,1) 1])) + (1-GP_params.nugget)*Full_Kcor;
                variational_params(i_cell).shapes.Full_Kcor=Full_Kcor;
                
                    sigma_mat=sqrt(variational_params(i_cell).shapes.vi_sigma)*ones(1,length(variational_params(i_cell).shapes.vi_sigma));
                    Full_Kcov= sigma_mat.*Full_Kcor.*sigma_mat';
                    variational_params(i_cell).shapes.Sigma_inv=inv(Full_Kcov);
                    variational_params(i_cell).shapes.Sigma_inv=( variational_params(i_cell).shapes.Sigma_inv+ variational_params(i_cell).shapes.Sigma_inv')/2;
                    variational_params(i_cell).shapes.Sigma=Full_Kcov;
%                     Dmat= diag( diag(variational_params(i_cell).shapes.Sigma_inv).*exp(-variational_params(i_cell).shapes.log_sigma) );
                    Dmat= diag( diag(variational_params(i_cell).shapes.Sigma_inv).*0);
                    variational_params(i_cell).shapes.Sigma_tilde_inv=variational_params(i_cell).shapes.Sigma_inv+Dmat;
                    variational_params(i_cell).shapes.Sigma_tilde=inv(variational_params(i_cell).shapes.Sigma_tilde_inv);
                    variational_params(i_cell).shapes.Sigma_tilde=(variational_params(i_cell).shapes.Sigma_tilde+variational_params(i_cell).shapes.Sigma_tilde')/2;
            end
        end
        
        prior_params=variational_params;
        % Assign the correct variance for prior distributions:
        for i_cell = 1:n_cell
            prior_params(i_cell).shapes.log_sigma=log(prior_params(i_cell).shapes.prior_sigma);
            sigma_mat=sqrt(prior_params(i_cell).shapes.prior_sigma)*ones(1,length(prior_params(i_cell).shapes.prior_sigma));
            Full_Kcov= sigma_mat.*prior_params(i_cell).shapes.Full_Kcor.*sigma_mat';
            prior_params(i_cell).shapes.Sigma=Full_Kcov;
            prior_params(i_cell).shapes.Sigma_tilde=Full_Kcov;
        end
    case 'single'
        % Single location for each neuron
        % One location corresponds to only one neuron!
        
%         boundary_params = 1; % allow for some errors
        for i_trial = 1:length(trials)
            stim_locs = trials(i_trial).locations;
            for i_loc = size(stim_locs,1)
                if ~isnan(stim_locs(i_loc,1))
                    trials(i_trial).cell_and_pos{i_loc}= zeros(0,2); %i_cell, i_pos
                    for i_cell = 1:n_cell
                        this_loc=neurons(i_cell).location;
                        rel_pos=stim_locs(i_loc,:)-this_loc;
%                         if sum(rel_pos.^2)<boundary_params
                            % Check for xy shape and z shape:
                            trials(i_trial).cell_and_pos{i_loc}=[trials(i_trial).cell_and_pos{i_loc};...
                                [i_cell 1]];
%                         end
                    end
                end
            end
        end
%         variational_params=rmfield(variational_params,{'shapes'});
        prior_params=variational_params;
end
%%
% GP_params=experiment_setup.prior_info.prior_parameters.GP_params