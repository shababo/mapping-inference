function stim_location_this_group = ...
    get_stim_locations(this_cell,group_ID,cell_params,design_params,template_cell)

% template_cell=experiment_setup.prior_info.template_cell;
cell_locations_neighbourhood = reshape([cell_params(:).location],3,[])';

stim_location_this_group=struct;

z_depth=mean(cell_locations_neighbourhood(:,3));
% The list of all related cells :
grid_coord = zeros(sum(design_params.candidate_grid_params.number)+1,3);
i_count=1;
grid_coord(i_count,:)=zeros(1,3);
for i_circle=1:length(design_params.candidate_grid_params.number)
    this_grid_number=design_params.candidate_grid_params.number(i_circle);
    for i_grid = 1:this_grid_number
        i_count=i_count+1;
        grid_coord(i_count,:)=...
            design_params.candidate_grid_params.radius(i_circle)*...
            [sin(2*pi*i_grid/this_grid_number) cos(2*pi*i_grid/this_grid_number) 0];
    end
end
grid_this_cell = grid_coord+ones(size(grid_coord,1),1)*cell_params(this_cell).location;


switch group_ID
    case 'undefined'
        
        grid_this_cell(:,3)= z_depth;
        
    case 'connected'
        
end
[effect_this_cell] = get_weights(cell_params, template_cell,grid_this_cell);
stim_location_this_group.grid=grid_this_cell;
stim_location_this_group.effect=effect_this_cell;

%--------------------------------------------------%
% Construct stim sets for the connected cells



