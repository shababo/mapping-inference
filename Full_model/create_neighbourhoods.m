function neighbourhoods = create_neighbourhoods(experiment_setup)
% initialize the neighbourhoods 

switch experiment_setup.experiment_type
    case 'simulation'
        
        cell_locations=reshape([experiment_setup.neurons(:).location]',[],length(experiment_setup.neurons(1).location));
        z_min=min(cell_locations(:,3));
        z_max=max(cell_locations(:,3));
        z_slide_width=experiment_setup.neighbourhood_params.height;
        z_borders=(z_min):z_slide_width:(z_max);
        number_of_neighbourhoods = length(z_borders);
        
        number_of_cells=size(cell_locations,1);
        cell_group_idx = zeros(number_of_cells,1);
        for i_cell = 1:size(cell_locations,1)
            cell_group_idx(i_cell)= sum(cell_locations(i_cell,3)>z_borders);
        end
        
        % Initialize the neighbourhoods
        % copying the neuron info from experiment setup 
        neighbourhoods=struct([]);
        for i_neighbourhood = 1:number_of_neighbourhoods
           %neighbourhoods(i_neighbourhood)=struct;
           neighbourhoods(i_neighbourhood).neighbourhood_ID=i_neighbourhood;
           %neighbourhoods(i_neighbourhood).neurons=struct;
           cell_list_this_neighbourhood=find(cell_group_idx==i_neighbourhood);
           for i_cell = 1:length(cell_list_this_neighbourhood)
               cell_ID=cell_list_this_neighbourhood(i_cell);
                neighbourhoods(i_neighbourhood).neurons(i_cell)=experiment_setup.neurons(cell_ID);
            end
           neighbourhoods(i_neighbourhood).computing_time=struct;
        end
        
        
               
    
        % Calculate the candidate grid for each cell in each neigbourhood 
        group_names = fieldnames(experiment_setup.groups);
        for i_neighbourhood = 1:number_of_neighbourhoods
            cell_locations_neighbourhood = reshape([neighbourhoods(i_neighbourhood).neurons(:).location]',[],3);
            for i_cell = 1:length(neighbourhoods(i_neighbourhood).neurons)
                neighbourhoods(i_neighbourhood).neurons(i_cell).stim_locations=struct([]);
                 neighbourhoods(i_neighbourhood).neurons(i_cell).cell_ID=cell_ID;
                neighbourhoods(i_neighbourhood).neurons(i_cell).group_ID='undefined'; % all are initialized as undefined
                neighbourhoods(i_neighbourhood).neurons(i_cell).primary_indicator=true;
              
                for i_group = 1:length(group_names)
                    neighbourhoods(i_neighbourhood).neurons(i_cell).stim_locations
                    
                    design_params=experiment_setup.groups.(group_names{i_group}).design_func_params;
                    cell_params=neighbourhoods(i_neighbourhood).neurons;
                    =[];...
                        get_stim_locations(cell_params,design_params);
                    % define grid, effect, and inner product 
                end
            end
        end
        
        % Initialize PR and gain parameters 
        for i_neighbourhood = 1:number_of_neighbourhoods
            for i_cell = 1:length(neighbourhoods(i_neighbourhood).neurons)
                neighbourhoods(i_neighbourhood).neurons(i_cell).PR_params=struct([]);
                neighbourhoods(i_neighbourhood).neurons(i_cell).PR_params(1).type='SpikedLogitNormal';
                neighbourhoods(i_neighbourhood).neurons(i_cell).PR_params(1).pi=0;
                neighbourhoods(i_neighbourhood).neurons(i_cell).PR_params(1).alpha=0;
                neighbourhoods(i_neighbourhood).neurons(i_cell).PR_params(1).beta=0;
                neighbourhoods(i_neighbourhood).neurons(i_cell).PR_params(1).mean=0;
                neighbourhoods(i_neighbourhood).neurons(i_cell).PR_params(1).variance=0;
                neighbourhoods(i_neighbourhood).neurons(i_cell).PR_params(1).upper_quantile=0;
                neighbourhoods(i_neighbourhood).neurons(i_cell).PR_params(1).lower_quantile=0;
                neighbourhoods(i_neighbourhood).neurons(i_cell).PR_params(1).nonzero_prob=0;
                
                %neighbourhoods(i_neighbourhood).neurons(i_cell).gain=struct([]);
            end
        end
        
    
    case 'experiment'
       
% REDO WITH NEW STRUCTURES!!!!! NOT DONE YET
% MAKES SURE TO PUT ALL CELLS IN INIT ANALYSIS GROUP (PROB. UNDEFINED)
if ~isfield(params.exp,'foe_bounds')
    params.exp.foe_bounds = [-148 148; -148 148];
end

experiment_setup.neurons 

neighbourhoods.cell_locations = cell_locations;
neighbourhoods.z_locs = z_locs;
neighbourhoods.n_cell = size(cell_locations,1);

% break up into z groups
z_slice_width = params.exp.z_width;
z_borders = [z_locs(1)-floor(z_slice_width/2); z_locs+floor(z_slice_width/2)];

neighbourhoods.n_planes = length(z_locs);

cell_group_idx = zeros(neighbourhoods.n_cell,1);
for i_cell = 1:size(cell_locations,1)
    cell_group_idx(i_cell)= sum(cell_locations(i_cell,3)>z_borders);
end
neighbourhoods.cell_group_list = cell(neighbourhoods.n_planes,1);
for i_plane = 1:neighbourhoods.n_planes
    neighbourhoods.cell_group_list{i_plane} = find(cell_group_idx==i_plane);
end

assignin('base','cell_group_list',neighbourhoods.cell_group_list)
assignin('base','cell_locations',cell_locations)

pi_target_selected = cell(neighbourhoods.n_planes,1);
inner_normalized_products = cell(neighbourhoods.n_planes,1);
target_locations_selected = cell(neighbourhoods.n_planes,1);
power_selected = cell(neighbourhoods.n_planes,1);
target_locations_all = cell(neighbourhoods.n_planes,1);
cell_neighbours = cell(neighbourhoods.n_planes,1);
target_locations_nuclei = cell(neighbourhoods.n_planes,1);
power_nuclei = cell(neighbourhoods.n_planes,1);
pi_target_nuclei = cell(neighbourhoods.n_planes,1);
loc_to_cell_nuclei = cell(neighbourhoods.n_planes,1);
loc_to_cell = cell(neighbourhoods.n_planes,1);

for i = 1:neighbourhoods.n_planes

    n_cell_this_plane = length(neighbourhoods.cell_group_list{i});
    target_cell_list.primary = neighbourhoods.cell_group_list{i};
    target_cell_list.secondary = [];
    [pi_target_selected{i}, inner_normalized_products{i},targ_locs_tmp,...
    targ_loc_nuc_tmp, pi_target_nuclei{i}, loc_to_cell_nuclei{i}] = ...
        get_stim_locations(...
        target_cell_list,cell_locations,...
        params.r1,params.r2,params.num_per_grid,params.num_per_grid_dense,params.template_cell.shape_template,...
        params.design.stim_loc_type,z_locs(i),params.exp.arbitrary_z);
    
    targ_locs_tmp(targ_locs_tmp(:,1) < params.exp.foe_bounds(1,1),1) = params.exp.foe_bounds(1,1);
    targ_locs_tmp(targ_locs_tmp(:,1) > params.exp.foe_bounds(1,2),1) = params.exp.foe_bounds(1,2);
    targ_locs_tmp(targ_locs_tmp(:,2) < params.exp.foe_bounds(2,1),2) = params.exp.foe_bounds(2,1);
    targ_locs_tmp(targ_locs_tmp(:,2) > params.exp.foe_bounds(2,2),2) = params.exp.foe_bounds(2,2);
    
    targ_loc_nuc_tmp(targ_loc_nuc_tmp(:,1) < params.exp.foe_bounds(1,1),1) = params.exp.foe_bounds(1,1);
    targ_loc_nuc_tmp(targ_loc_nuc_tmp(:,1) > params.exp.foe_bounds(1,2),1) = params.exp.foe_bounds(1,2);
    targ_loc_nuc_tmp(targ_loc_nuc_tmp(:,2) < params.exp.foe_bounds(2,1),2) = params.exp.foe_bounds(2,1);
    targ_loc_nuc_tmp(targ_loc_nuc_tmp(:,2) > params.exp.foe_bounds(2,2),2) = params.exp.foe_bounds(2,2);
    
    target_locations_selected{i} = targ_locs_tmp;
    target_locations_nuclei{i} = targ_loc_nuc_tmp;
    
    loc_to_cell{i} = zeros(size(pi_target_selected,2),1);
    for i_cell = 1:length(target_cell_list.primary)
        loc_to_cell{i}( (i_cell-1)*(2*params.num_per_grid+1)+ (1:(2*params.num_per_grid+1)))=i_cell;     
    end  
end




neighbourhoods.pi_target_selected = pi_target_selected;
neighbourhoods.inner_normalized_products = inner_normalized_products;
neighbourhoods.target_locations_selected = target_locations_selected;
neighbourhoods.target_locations_nuclei = target_locations_nuclei;
neighbourhoods.pi_target_nuclei = pi_target_nuclei;
neighbourhoods.loc_to_cell_nuclei = loc_to_cell_nuclei;
neighbourhoods.loc_to_cell = loc_to_cell;
    case 'reproduction'
        
end
