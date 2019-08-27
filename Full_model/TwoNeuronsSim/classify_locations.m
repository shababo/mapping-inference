function [new_shape_params,new_trials]=classify_locations(current_params,neurons,new_trials,prior_info)
%% Output:
% new_shape_params: contain basic information of the new locations relative
% to each neuron
% new_trials: modify the indicate whether the locations are new

%% Detect new locations:

GP_params=prior_info.GP_params;
clear('new_shape_params');
min_dist = GP_params.min_dist;
n_cell=length(neurons);
for i_cell = 1:n_cell
    new_shape_params(i_cell)=neurons(i_cell).params(1);
end
%% Calculate the parameters in the conditional distribution (correlation matrix for the shapes)


boundary_params = GP_params.GP_boundary;
type=GP_params.type;
if strcmp(type, 'xy_square')
    n_axis_model =2;
    axis_names = {'xy', 'z'};
    for i_trial = 1:length(new_trials)
        stim_locs = new_trials(i_trial).locations;
        for i_loc = size(stim_locs,1)
            if ~isnan(stim_locs(i_loc,1))
                new_trials(i_trial).cell_and_pos{i_loc}= zeros(0,3); %i_cell, i_pos
                new_trials(i_trial).cell_and_status{i_loc}= zeros(0,3); %i_cell, i_pos
                for i_cell = 1:n_cell
                    this_loc=neurons(i_cell).location;
                    rel_pos=stim_locs(i_loc,:)-this_loc;
                    if check_in_boundary(rel_pos,boundary_params)
                        % Check for xy shape and z shape:
                        ia=zeros(1,2);
                        status_flag=zeros(1,2); % 0: a new position; 1: same as in the fitted locations;
                        
                        for i_axis = 1:n_axis_model
                            if i_axis == 1
                                this_rel_pos=rel_pos(1:2);
                            else
                                this_rel_pos=rel_pos(3);
                            end
                            % Instead of exact matching, find if the new
                            % stim loc is within 2 microns from existing
                            % ones (fitted)
                            existing_loc_dim=size(current_params(i_cell).(axis_names{i_axis}).locations);
                            sq_dist=(current_params(i_cell).(axis_names{i_axis}).locations- ones(existing_loc_dim(1),1)*this_rel_pos).^2;
                            if existing_loc_dim(2) > 1
                                sq_dist = sum(sq_dist,2);
                            end
                            if (existing_loc_dim(1) == 0) | min(sq_dist.^(1/2))> min_dist
                                tmp = [];
                            else
                                [~,tmp]=min(sq_dist.^(1/2));status_flag(i_axis)=1;
                            end
                            %[C,tmp,~] = intersect(variational_params(i_cell).(axis_names{i_axis}).locations,this_rel_pos,'rows');
                            if status_flag(i_axis)==0
                                new_loc_dim=size(new_shape_params(i_cell).(axis_names{i_axis}).locations);
                                if (new_loc_dim(1) > 0)
                                    sq_dist=(new_shape_params(i_cell).(axis_names{i_axis}).locations- ones(new_loc_dim(1),1)*this_rel_pos).^2;
                                
                                if new_loc_dim(2) > 1
                                    
                                    sq_dist = sum(sq_dist,2);
                                end
                                end 
                                if (new_loc_dim(1) == 0) | min(sq_dist.^(1/2))> min_dist
                                    tmp = [];
                                else
                                    [~,tmp]=min(sq_dist.^(1/2));status_flag(i_axis)=0;
                                end
                            end
                            if isempty(tmp) % this is a brand new loc
                                new_shape_params(i_cell).(axis_names{i_axis}).locations=...
                                    [new_shape_params(i_cell).(axis_names{i_axis}).locations; this_rel_pos];
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
                                new_shape_params(i_cell).(axis_names{i_axis}).bounds.low = [new_shape_params(i_cell).(axis_names{i_axis}).bounds.low; lower_bound];
                                new_shape_params(i_cell).(axis_names{i_axis}).bounds.up = [new_shape_params(i_cell).(axis_names{i_axis}).bounds.up; upper_bound];
                                % logit transform:
                                switch  new_shape_params(i_cell).(axis_names{i_axis}).dist
                                    case 'logit-normal'
                                        mean_logit=log( (mean_3d-lower_bound)/(upper_bound -mean_3d)); % this is actually 0
                                        new_shape_params(i_cell).(axis_names{i_axis}).mean=[new_shape_params(i_cell).(axis_names{i_axis}).mean; mean_logit];
                                        new_shape_params(i_cell).(axis_names{i_axis}).log_sigma=[new_shape_params(i_cell).(axis_names{i_axis}).log_sigma; 0];% variance is no longer the original one!
                                    case 'mvn'
                                        new_shape_params(i_cell).(axis_names{i_axis}).prior_sigma=[new_shape_params(i_cell).(axis_names{i_axis}).prior_sigma; ...
                                            var_3d];
                                        %                                          variational_params(i_cell).shapes.mean=[variational_params(i_cell).(axis_names{i_axis}).mean; 0];
                                        new_shape_params(i_cell).(axis_names{i_axis}).mean=[new_shape_params(i_cell).(axis_names{i_axis}).mean; ...
                                            log( (mean_3d-lower_bound)/(upper_bound-mean_3d))];
                                        %variational_params(i_cell).(axis_names{i_axis}).mean=[variational_params(i_cell).(axis_names{i_axis}).mean; ...
                                        %  mean_3d];
                                        new_shape_params(i_cell).(axis_names{i_axis}).log_sigma=log(new_shape_params(i_cell).(axis_names{i_axis}).prior_sigma);
                                        %variational_params(i_cell).(axis_names{i_axis}).log_sigma(:)=2;
                                end
                                ia(i_axis)=length(new_shape_params(i_cell).(axis_names{i_axis}).mean);
                            else
                                ia(i_axis)=tmp;
                            end
                            
                        end
                        new_trials(i_trial).cell_and_pos{i_loc}=[new_trials(i_trial).cell_and_pos{i_loc};...
                            i_cell ia];
                        new_trials(i_trial).cell_and_status{i_loc}=[new_trials(i_trial).cell_and_status{i_loc};...
                            i_cell status_flag];
                        
                    end
                end
            end
        end
    end
end

