function stim_location_this_group = ...
    get_stim_locations(this_cell,group_ID,cell_params,design_params,varargin)

% this_cell=i_cell;
% group_ID=group_names{i_group};
% template_cell=experiment_setup.prior_info.template_cell;
% foe_bounds=experiment_setup.exp.foe_bounds;
% z_depth = mean(cell_locations_neighbourhood(:,3));
if length(varargin) > 1 && ~isempty(varargin{2})
    foe_bounds = varargin{2};
else
    foe_bounds = [];
end

% template_cell=experiment_setup.prior_info.template_cell;
cell_locations_neighbourhood = reshape([cell_params(:).location],3,[])';
if ~isempty(varargin) && ~isempty(varargin{1})
    z_depth = varargin{1};
else
    z_depth = mean(cell_locations_neighbourhood(:,3));
end
stim_location_this_group=struct;


% The list of all related cells :
grid_coord = zeros(sum(design_params.candidate_grid_params.number)+1,3);
i_count=1;
grid_coord(i_count,:)=zeros(1,3);
for i_circle=1:length(design_params.candidate_grid_params.radius)
    this_grid_number=design_params.candidate_grid_params.number(i_circle);
    for i_grid = 1:this_grid_number
        i_count=i_count+1;
        grid_coord(i_count,:)=...
            design_params.candidate_grid_params.radius(i_circle)*...
            [sin(2*pi*i_grid/this_grid_number) cos(2*pi*i_grid/this_grid_number) 0];
    end
end
grid_this_cell = grid_coord+ones(size(grid_coord,1),1);


switch group_ID
    case 'undefined'
        
        grid_this_cell(:,3)= z_depth;
        
    case 'connected'
        
    case 'alive'
        
end

% cut based on slm range
if ~isempty(foe_bounds)
    for i_dim = 1:size(foe_bounds,1)
        grid_this_cell(grid_this_cell(:,i_dim) < foe_bounds(i_dim,1),i_dim) = foe_bounds(i_dim,1);
        grid_this_cell(grid_this_cell(:,i_dim) > foe_bounds(i_dim,2),i_dim) = foe_bounds(i_dim,2);
    end
end


% [effect_this_cell] = get_weights(cell_params, template_cell,grid_this_cell);
stim_location_this_group.grid=grid_this_cell;
% stim_location_this_group.effect=effect_this_cell;

%--------------------------------------------------%
% Construct stim sets for the connected cells



