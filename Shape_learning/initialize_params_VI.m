function [variational_params, prior_params,trials]=initialize_params_VI(neurons,trials,prior_info)
%%
n_cell=length(neurons);
clear('variational_params')

GP_params=prior_info.GP_params;
type=GP_params.type;
min_dist = 2; % minimun distance to distinguish two stim spot 
%variational_params(n_cell)=struct;
for i_cell = 1:n_cell
    variational_params(i_cell)=neurons(i_cell).params(end);
    % Check if there are new locations in this batch
end

boundary_params = prior_info.GP_params.GP_boundary;


if strcmp(type, 'xy_square')
    n_axis_model =2;
    axis_names = {'xy', 'z'};
    
    for i_trial = 1:length(trials)
        stim_locs = trials(i_trial).locations;
        for i_loc = size(stim_locs,1)
            if ~isnan(stim_locs(i_loc,1))
                trials(i_trial).cell_and_pos{i_loc}= zeros(0,3); %i_cell, i_pos
                for i_cell = 1:n_cell
                    this_loc=neurons(i_cell).location;
                    rel_pos=stim_locs(i_loc,:)-this_loc;
                    if check_in_boundary(rel_pos,boundary_params)
                        % Check for xy shape and z shape:
                        ia=zeros(1,2);
                        
                        for i_axis = 1:n_axis_model
                            if i_axis == 1
                                this_rel_pos=rel_pos(1:2);
                            else
                                this_rel_pos=rel_pos(3);
                            end
                            
                            % Instead of exact matching, find if the new
                            % stim loc is within 2 microns from existing
                            % ones
                            existing_loc_dim=size(variational_params(i_cell).(axis_names{i_axis}).locations);
                            sq_dist=(variational_params(i_cell).(axis_names{i_axis}).locations- ones(existing_loc_dim(1),1)*this_rel_pos).^2;
                            if existing_loc_dim(2) > 1
                                sq_dist = sum(sq_dist,2);
                            end
                            if (existing_loc_dim(1) == 0) | min(sq_dist.^(1/2))> min_dist
                                tmp = []; 
                            else
                                
                                 [~,tmp]=min(sq_dist.^(1/2));
                            end
                            %[C,tmp,~] = intersect(variational_params(i_cell).(axis_names{i_axis}).locations,this_rel_pos,'rows');
                            
                            if isempty(tmp) % this is a new location:
                                
                                variational_params(i_cell).(axis_names{i_axis}).locations=...
                                    [variational_params(i_cell).(axis_names{i_axis}).locations; this_rel_pos];
                                
                                [interpolated_shape]=interpolate_3D(rel_pos,GP_params,type);
                                if i_axis == 1
                                    mean_3d=interpolated_shape.mean_xy;
                                    var_3d=interpolated_shape.var_xy;
                                else
                                    
                                    mean_3d=interpolated_shape.mean_z;
                                    var_3d=interpolated_shape.var_z;
                                end
                                % Adding a small variance to allow the shape to
                                % learn from the data set when the prior shape
                                % variance is too small.
                                if (var_3d <  prior_info.GP_params.GP_minimal_variance)
                                    var_3d =   prior_info.GP_params.GP_minimal_variance;
                                end
                                if isfield(prior_info, 'GP_added_variance')
                                    var_3d = var_3d+prior_info.GP_params.GP_added_variance;
                                end
                                lower_bound =max(0, mean_3d-2*sqrt(var_3d));upper_bound =min(1, mean_3d+2*sqrt(var_3d));
%                                 lower_bound =0;upper_bound =1;
                                variational_params(i_cell).(axis_names{i_axis}).bounds.low = [variational_params(i_cell).(axis_names{i_axis}).bounds.low; lower_bound];
                                variational_params(i_cell).(axis_names{i_axis}).bounds.up = [variational_params(i_cell).(axis_names{i_axis}).bounds.up; upper_bound];
                                % logit transform:
                                switch  variational_params(i_cell).(axis_names{i_axis}).dist
                                    case 'logit-normal'
                                        mean_logit=log( (mean_3d-lower_bound)/(upper_bound -mean_3d)); % this is actually 0
                                        variational_params(i_cell).(axis_names{i_axis}).mean=[variational_params(i_cell).(axis_names{i_axis}).mean; mean_logit];
                                        variational_params(i_cell).(axis_names{i_axis}).log_sigma=[variational_params(i_cell).(axis_names{i_axis}).log_sigma; 0];% variance is no longer the original one!
                                    case 'mvn'
                                        variational_params(i_cell).(axis_names{i_axis}).prior_sigma=[variational_params(i_cell).(axis_names{i_axis}).prior_sigma; ...
                                           var_3d];
%                                          variational_params(i_cell).shapes.mean=[variational_params(i_cell).(axis_names{i_axis}).mean; 0];
                                        variational_params(i_cell).(axis_names{i_axis}).mean=[variational_params(i_cell).(axis_names{i_axis}).mean; ...
                                           log( (mean_3d-lower_bound)/(upper_bound-mean_3d))];
                                         %variational_params(i_cell).(axis_names{i_axis}).mean=[variational_params(i_cell).(axis_names{i_axis}).mean; ...
                                          %  mean_3d];
                                        variational_params(i_cell).(axis_names{i_axis}).log_sigma=log(variational_params(i_cell).(axis_names{i_axis}).prior_sigma);
                                        variational_params(i_cell).(axis_names{i_axis}).log_sigma(:)=2;
                                end
                                ia(i_axis)=length(variational_params(i_cell).(axis_names{i_axis}).mean);
                            else
                                ia(i_axis)=tmp;
                            end
                            
                        end
                        trials(i_trial).cell_and_pos{i_loc}=[trials(i_trial).cell_and_pos{i_loc};...
                            i_cell ia];
                    end
                end
            end
        end
    end
    % Calculate the joint distribution (correlation matrix for the shapes)
    if strcmp(variational_params(1).(axis_names{i_axis}).dist,'mvn')
        for i_cell = 1:n_cell
            for i_axis = 1:n_axis_model
                
                X=variational_params(i_cell).(axis_names{i_axis}).locations;
                
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
                    tmp_Kcor=get_kernel_cor(X,X,tau);
                    tmp_Kcor=(tmp_Kcor+tmp_Kcor')/2;
                    Full_Kcor=tmp_Kcor;
                    
                end
                sigma_mat=variational_params(i_cell).(axis_names{i_axis}).prior_sigma*ones(1,length(variational_params(i_cell).(axis_names{i_axis}).prior_sigma));
                Full_Kcov= sigma_mat.*Full_Kcor.*sigma_mat';
                
                variational_params(i_cell).(axis_names{i_axis}).Sigma_inv=inv(Full_Kcov);
                variational_params(i_cell).(axis_names{i_axis}).Sigma_inv=( variational_params(i_cell).(axis_names{i_axis}).Sigma_inv+ variational_params(i_cell).(axis_names{i_axis}).Sigma_inv')/2;
                variational_params(i_cell).(axis_names{i_axis}).Sigma=Full_Kcov;
                Dmat= diag( diag(variational_params(i_cell).(axis_names{i_axis}).Sigma_inv).*exp(-variational_params(i_cell).(axis_names{i_axis}).log_sigma) );
                
                variational_params(i_cell).(axis_names{i_axis}).Sigma_tilde_inv=variational_params(i_cell).(axis_names{i_axis}).Sigma_inv+Dmat;
                variational_params(i_cell).(axis_names{i_axis}).Sigma_tilde=inv(variational_params(i_cell).(axis_names{i_axis}).Sigma_inv);
                variational_params(i_cell).(axis_names{i_axis}).Sigma_tilde=(variational_params(i_cell).(axis_names{i_axis}).Sigma_tilde+variational_params(i_cell).(axis_names{i_axis}).Sigma_tilde')/2;
            end
        end
        
    end
    prior_params=variational_params;
else
    
    
    for i_trial = 1:length(trials)
        stim_locs = trials(i_trial).locations;
        for i_loc = size(stim_locs,1)
            if ~isnan(stim_locs(i_loc,1))
                trials(i_trial).cell_and_pos{i_loc}= zeros(0,2); %i_cell, i_pos
                
                for i_cell = 1:n_cell
                    this_loc=neurons(i_cell).location;
                    rel_pos=stim_locs(i_loc,:)-this_loc;
                    if check_in_boundary(rel_pos,boundary_params)
                        [C,ia,ib] = intersect(variational_params(i_cell).shapes.locations,rel_pos,'rows');
                        if isempty(C) % this is a new location:
                            variational_params(i_cell).shapes.locations=...
                                [variational_params(i_cell).shapes.locations; rel_pos];
                            
                            [interpolated_shape]=interpolate_3D(rel_pos,GP_params,type);
                            mean_3d=interpolated_shape.mean_3d;
                            var_3d=interpolated_shape.var_3d;
                            % Adding a small variance to allow the shape to
                            % learn from the data set when the prior shape
                            % variance is too small.
                            if (var_3d <   prior_info.prior_parameters.GP_minimal_variance)
                                var_3d =   prior_info.prior_parameters.GP_minimal_variance;
                            end
                            if isfield(prior_info, 'GP_added_variance')
                                var_3d = var_3d+prior_info.GP_added_variance;
                            end
                            %                         lower_bound =max(0, mean_3d-4*sqrt(var_3d));upper_bound =min(1, mean_3d+4*sqrt(var_3d));
                            lower_bound =0;upper_bound =1;
                            variational_params(i_cell).shapes.bounds.low = [variational_params(i_cell).shapes.bounds.low; lower_bound];
                            
                            variational_params(i_cell).shapes.bounds.up = [variational_params(i_cell).shapes.bounds.up; upper_bound];
                            
                            % logit transform:
                            switch  variational_params(i_cell).shapes.dist
                                case 'logit-normal'
                                    mean_logit=log( (mean_3d-lower_bound)/(upper_bound -mean_3d)); % this is actually 0
                                    variational_params(i_cell).shapes.mean=[variational_params(i_cell).shapes.mean; mean_logit];
                                    variational_params(i_cell).shapes.log_sigma=[variational_params(i_cell).shapes.log_sigma; 0];% variance is no longer the original one!
                                case 'mvn'
                                    variational_params(i_cell).shapes.prior_sigma=[variational_params(i_cell).shapes.prior_sigma; sqrt(var_3d)];
                                    %                                 variational_params(i_cell).shapes.mean=[variational_params(i_cell).shapes.mean; 0];
                                    variational_params(i_cell).shapes.mean=[variational_params(i_cell).shapes.mean; log(mean_3d/(1-mean_3d))];
                                    
                                    variational_params(i_cell).shapes.log_sigma=log(variational_params(i_cell).shapes.prior_sigma);
                                    variational_params(i_cell).shapes.log_sigma(:)=2;
                            end
                            trials(i_trial).cell_and_pos{i_loc}=[trials(i_trial).cell_and_pos{i_loc};...
                                i_cell length(variational_params(i_cell).shapes.mean)];
                        else
                            trials(i_trial).cell_and_pos{i_loc}=[trials(i_trial).cell_and_pos{i_loc};...
                                i_cell ia];
                        end
                    end
                    
                end
            end
        end
    end
    
    % Calculate the joint distribution (correlation matrix for the shapes)
    if strcmp(variational_params(1).shapes.dist,'mvn')
        for i_cell = 1:n_cell
            locs=variational_params(i_cell).shapes.locations;
            
            axis_list= fieldnames(GP_params);
            
            tmp_i=strcmp(axis_list,'type');
%             axis_list=axis_list(find(~tmp_i));
            axis_list={'x', 'y', 'z'};
            
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
            sigma_mat=variational_params(i_cell).shapes.prior_sigma*ones(1,length(variational_params(i_cell).shapes.prior_sigma));
            Full_Kcov= sigma_mat.*Full_Kcor.*sigma_mat';
            
            variational_params(i_cell).shapes.Sigma_inv=inv(Full_Kcov);
            variational_params(i_cell).shapes.Sigma_inv=( variational_params(i_cell).shapes.Sigma_inv+ variational_params(i_cell).shapes.Sigma_inv')/2;
            variational_params(i_cell).shapes.Sigma=Full_Kcov;
            Dmat= diag( diag(variational_params(i_cell).shapes.Sigma_inv).*exp(-variational_params(i_cell).shapes.log_sigma) );
            
            variational_params(i_cell).shapes.Sigma_tilde_inv=variational_params(i_cell).shapes.Sigma_inv+Dmat;
            variational_params(i_cell).shapes.Sigma_tilde=inv(variational_params(i_cell).shapes.Sigma_inv);
            variational_params(i_cell).shapes.Sigma_tilde=(variational_params(i_cell).shapes.Sigma_tilde+variational_params(i_cell).shapes.Sigma_tilde')/2;
        end
    end
    prior_params=variational_params;
end
%%

% GP_params=experiment_setup.prior_info.prior_parameters.GP_params