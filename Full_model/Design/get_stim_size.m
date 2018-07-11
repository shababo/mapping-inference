function [stimuli_size] = get_stim_size(trials,neighbourhood,experiment_setup)
%% Debug:
% neighbourhood = neighbourhoods(1);
% trials = experiment_query.undefined.trials;
%% Dimensions
number_of_trials=length(trials);
number_of_cells =length(neighbourhood.neurons);
number_of_spots = size(trials(1).locations,1);
simulation_params=experiment_setup.sim;

boundary_params = experiment_setup.prior_info.prior_parameters.boundary_params;
%%
stimuli_size=zeros(number_of_trials,number_of_cells);
cell_ID_list=[neighbourhood.neurons(:).cell_ID];

for l = 1:number_of_trials
    for m = 1:number_of_spots
        %         this_loc_ID=trials(l).location_IDs(m);
        if m>size(trials(l).locations,1)
        else
            this_loc=trials(l).locations(m,:);
            
            this_loc_power =trials(l).power_levels(m);
            %             cell_ID=trials(l).cell_IDs(m);
            for i_cell = 1:number_of_cells
                
                rel_loc = this_loc - neighbourhood.neurons(i_cell).truth.location;
                
                if (neighbourhood.neurons(i_cell).truth.PR & check_in_boundary(rel_loc, boundary_params))>0
                    %             fprintf([num2str(l) '\n'])
                    this_size = griddata(simulation_params.mesh_grid(:,1),simulation_params.mesh_grid(:,2),simulation_params.mesh_grid(:,3),...
                        neighbourhood.neurons(i_cell).truth.shape,rel_loc(1),rel_loc(2),rel_loc(3));
                    if isnan(this_size)
                        this_size = 0;
                    end
                    stimuli_size(l,i_cell) = stimuli_size(l,i_cell)+ this_loc_power*this_size;
                end
            end
        end
    end
end


