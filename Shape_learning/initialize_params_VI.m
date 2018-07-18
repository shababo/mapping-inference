function [variational_params, prior_params,trials]=initialize_params_VI(neurons,trials,prior_info)
%%
n_cell=length(neurons);
clear('variational_params')
type='square';
GP_params=prior_info.prior_parameters.GP_params;
%variational_params(n_cell)=struct;
for i_cell = 1:n_cell
    variational_params(i_cell)=neurons(i_cell).params(end);
    % Check if there are new locations in this batch
    
end

boundary_params = prior_info.prior_parameters.boundary_params;
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
                        [mean_3d, var_3d]=interpolate_3D(rel_pos,GP_params,type);
                        lower_bound =max(0, mean_3d-2*sqrt(var_3d));upper_bound =min(1, mean_3d+2*sqrt(var_3d));
                        variational_params(i_cell).shapes.bounds.low = [variational_params(i_cell).shapes.bounds.low; lower_bound];
                        variational_params(i_cell).shapes.bounds.up = [variational_params(i_cell).shapes.bounds.up; upper_bound];
                        
                        % logit transform:
                        mean_logit=log( (mean_3d-lower_bound)/(upper_bound -mean_3d)); % this is actually 0
                        variational_params(i_cell).shapes.mean=[variational_params(i_cell).shapes.mean; mean_logit];
                        variational_params(i_cell).shapes.log_sigma=[variational_params(i_cell).shapes.log_sigma; 0]; % variance is no longer the original one!
%                         variational_params(i_cell).shapes.mean=[variational_params(i_cell).shapes.mean; mean_3d];
%                         variational_params(i_cell).shapes.log_sigma=[variational_params(i_cell).shapes.log_sigma; log(sqrt(var_3d))]; % for marginalization
                        
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


prior_params=variational_params;
