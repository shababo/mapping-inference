 function [variational_params, prior_params]=initialize_params_VI(neurons,trials,prior_info)
    n_cell=length(neurons);
    clear('variational_params')
    type='square';
    GP_params=prior_info.prior_parameters.GP_params;
    %variational_params(n_cell)=struct;
    for i_cell = 1:n_cell
        variational_params(i_cell)=neurons(i_cell).params(end);
        % Check if there are new locations in this batch
        this_loc=neurons(i_cell).location;
        for i_trial = 1:length(trials)
            stim_locs = trials(i_trial).locations;
            for i_loc = size(stim_locs,1)
               if ~isnan(stim_locs(i_loc,1))
                  rel_loc=stim_locs(i_loc,:)-this_loc; 
                   [C,ia,ib] = intersect(variational_params(i_cell).shapes.locations,rel_loc,'rows');
                   if isempty(C) % this is a new location:
                       variational_params(i_cell).shapes.locations=...
                           [variational_params(i_cell).shapes.locations; rel_loc];
                       [mean_3d, var_3d]=interpolate_3D(rel_loc,GP_params,type);
                       variational_params(i_cell).shapes.mean=[variational_params(i_cell).shapes.mean; mean_3d];
                       variational_params(i_cell).shapes.log_sigma=[variational_params(i_cell).shapes.log_sigma; var_3d];
                   end
               end
            end
        end 
    end
    
    
    prior_params=variational_params;
    