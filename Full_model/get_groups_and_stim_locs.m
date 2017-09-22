function data = get_groups_and_stim_locs(cell_locations,params,z_locs)

if ~isfield(params.exp,'foe_bounds')
    params.exp.foe_bounds = [-148 148; -148 148];
end

data.cell_locations = cell_locations;
data.z_locs = z_locs;
data.n_cell = size(cell_locations,1);

% break up into z groups
z_slice_width = params.exp.z_width;
z_borders = [z_locs(1)-floor(z_slice_width/2) z_locs+floor(z_slice_width/2)];

data.n_planes = length(z_locs);

cell_group_idx = zeros(data.n_cell,1);
for i_cell = 1:size(cell_locations,1)
    cell_group_idx(i_cell)= sum(cell_locations(i_cell,3)>z_borders);
end
data.cell_group_list = cell(data.n_planes,1);
for i_plane = 1:data.n_planes
    data.cell_group_list{i_plane} = find(cell_group_idx==i_plane);
end

r1=5;r2=10;r3=15;num_per_grid=12;

num_per_grid_dense = 16;
stim_threshold = params.eff_stim_threshold/params.template_cell.gain_template;

pi_target_selected = cell(data.n_planes,1);
inner_normalized_products = cell(data.n_planes,1);
target_locations_selected = cell(data.n_planes,1);
power_selected = cell(data.n_planes,1);
target_locations_all = cell(data.n_planes,1);
cell_neighbours = cell(data.n_planes,1);
target_locations_nuclei = cell(data.n_planes,1);
power_nuclei = cell(data.n_planes,1);
pi_target_nuclei = cell(data.n_planes,1);
loc_to_cell_nuclei = cell(data.n_planes,1);



stim_threshold = params.eff_stim_threshold/params.template_cell.gain_template;
for i = 1:data.n_planes

    n_cell_this_plane = length(data.cell_group_list{i});
    target_cell_list.primary = data.cell_group_list{i};
    target_cell_list.secondary = [];
    [pi_target_selected{i}, inner_normalized_products{i},targ_locs_tmp,power_selected{i},...
    target_locations_all{i},cell_neighbours{i},...
    targ_loc_nuc_tmp, power_nuclei{i},pi_target_nuclei{i}, loc_to_cell_nuclei{i}] = ...
        get_stim_locations(...
        target_cell_list,cell_locations,params.exp.user_power_level,...
        r1,r2,r3,num_per_grid,num_per_grid_dense,params.template_cell.shape_template,...
        params.stim_unique,params.template_cell.prob_trace,stim_threshold,params.design.stim_loc_type,z_locs(i),params.exp.arbitrary_z);
    
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
end




data.pi_target_selected = pi_target_selected;
data.inner_normalized_products = inner_normalized_products;
data.target_locations_selected = target_locations_selected;
data.power_selected = power_selected;
data.target_locations_all = target_locations_all;
data.cell_neighbours = cell_neighbours;
data.target_locations_nuclei = target_locations_nuclei;
data.power_nuclei = power_nuclei;
data.pi_target_nuclei = pi_target_nuclei;
data.loc_to_cell_nuclei = loc_to_cell_nuclei;