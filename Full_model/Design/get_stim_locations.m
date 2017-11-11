function [pi_target, inner_normalized_products,target_locations,loc_to_cell,...
    target_locations_nuclei,pi_target_nuclei, loc_to_cell_nuclei] = ...
    get_stim_locations(group_ID,cell_params,design_params,varargin)

cell_params
design_params

if ~isempty(varargin) && ~isempty(varargin{1})
    z_depth = varargin{1};
else
    cell_loc=reshape([cell_params(target_cell_list.primary).location],3,[])';
    z_depth = mean(cell_loc(:,3));
end

if length(varargin) > 1 && ~isempty(varargin{2})
    connected_arbitrary_z = varargin{2};
else
    connected_arbitrary_z = 0;
end

if length(varargin) > 1 && ~isempty(varargin{3})
    simulation_indicator = varargin{3};
else
    simulation_indicator = false;
end

% The list of all related cells :
cell_template = struct();
cell_template.shape= shape_template;

grid_coord = zeros(sum([grid_params(:).number])+1,3);
i_count=1;
grid_coord(i_count,:)=zeros(1,3);
for i_circle=1:length(grid_params)
    for i_grid = 1:grid_params(i_circle).number
        i_count=i_count+1;
        grid_coord(i_count,:)=...
            grid_params(i_circle).radius*[sin(2*pi*i_grid/grid_params(i_circle).number) cos(2*pi*i_grid/grid_params(i_circle).number) 0];
    end
end

% Calculate the stimulation locations 
target_locations = zeros(length(target_cell_list.primary)*size(grid_coord,1),3);
loc_to_cell=zeros(length(target_cell_list.primary)*size(grid_coord,1),1);
i_current=0;
for i_cell_index=1:length(target_cell_list.primary)
    i_cell= target_cell_list.primary(i_cell_index);
    grid_locs=cell_params(i_cell).location+grid_coord;
    target_locations(i_current+(1:size(grid_locs,1)),:) = grid_locs;
    loc_to_cell(i_current+(1:size(grid_locs,1)))=i_cell;
    i_current=i_current+size(grid_locs,1);
end
target_locations(:,3)= z_depth;

% [pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
%     cell_template,target_locations);
[pi_target, inner_normalized_products] = get_weights(related_cell_params, ...
    cell_template,target_locations,simulation_indicator);

%--------------------------------------------------%
% Construct stim sets for the connected cells


% Calculate the stimulation locations 
target_locations_nuclei = zeros(length(target_cell_list.primary)*size(grid_coord,1),3);
loc_to_cell_nuclei=zeros(length(target_cell_list.primary)*size(grid_coord,1),1);
i_current=0;
for i_cell_index=1:length(target_cell_list.primary)
    i_cell= target_cell_list.primary(i_cell_index);
    grid_locs=cell_params(i_cell).location+grid_coord;
    target_locations_nuclei(i_current+(1:size(grid_locs,1)),:) = grid_locs;
    loc_to_cell_nuclei(i_current+(1:size(grid_locs,1)))=i_cell;
    i_current=i_current+size(grid_locs,1);
end
if ~connected_arbitrary_z
    target_locations_nuclei(:,3) = z_depth;
end
[pi_target_nuclei, ~] = get_weights(related_cell_params, ...
    cell_template,target_locations_nuclei,simulation_indicator);

