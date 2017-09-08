function [pi_target_selected, inner_normalized_products,target_locations_selected,power_selected,...
    target_locations_all,cell_neighbours,...
 target_locations_3d_selected,power_3d_selected,pi_target_3d, inner_normalized_products_3d] = ...
    get_stim_locations(...
    cell_list,cell_locations,power_level,...
    r1,r2,r3,num_per_grid,num_per_grid_dense,num_z_grid,shape_template,...
    stim_unique,prob_trace,stim_threshold)

grid_jitters = zeros(num_per_grid,2);
for i_grid = 1:num_per_grid
   grid_jitters(i_grid,:)=[sin(2*pi*i_grid/num_per_grid) cos(2*pi*i_grid/num_per_grid)]; 
end
grid_jitters=[grid_jitters zeros(num_per_grid,1)];

% Calculate the stimulation locations 
target_locations = zeros(length(cell_list)*(2*num_per_grid+1),3);
for i_cell_index=1:length(cell_list)
    i_cell= cell_list(i_cell_index);
    nucleus_loc=cell_locations(i_cell,:);
    grid_locs=nucleus_loc;
    grid_locs=[grid_locs;bsxfun(@plus,grid_jitters*r1,nucleus_loc)];
    grid_locs=[grid_locs;bsxfun(@plus,grid_jitters*r2,nucleus_loc)];
    target_idx=(i_cell_index-1)*(2*num_per_grid+1) +(1: (2*num_per_grid+1));
    target_locations(target_idx,:) = grid_locs;
end
target_locations(:,3)= mean(cell_locations(cell_list,3));

%plot(target_locations{this_plane}(:,2),target_locations{this_plane}(:,1),'.')

cell_params.locations =  cell_locations(cell_list,:);
cell_params.shape_gain = ones(length(cell_list),1);
cell_template = struct();
cell_template.shape= shape_template;
% [pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
%     cell_template,target_locations);
[pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
    cell_template,target_locations);


%---- Select stim locations:

prob_cells_by_locations=zeros(size(pi_target,2),size(pi_target,1));
for l = 1:size(pi_target,2)
    for i_cell_index = 1:length(cell_list)
        prob_cells_by_locations(l,i_cell_index)= stim_to_prob(...
            pi_target(i_cell_index,l)*power_level(1),stim_unique,prob_trace);
    end
end

% Now, for each cell, find if there is a location that it is the only cell been stimulated 
unique_stim_locations = zeros(length(cell_list),1);
unique_stim_powers = zeros(length(cell_list),1);

num_cells_stimulated=sum(prob_cells_by_locations>0.01,2);
for i_cell_index = 1:length(cell_list)
    stim_for_this_cell=find( (num_cells_stimulated <2) & (prob_cells_by_locations(:,i_cell_index ) >0.8));
    if isempty(stim_for_this_cell)
        % do nothing
    else
        unique_stim_locations(i_cell_index)=min(stim_for_this_cell);
        unique_stim_powers(i_cell_index)=power_level(1);
    end
end
% Search for more stim locations at lower power level


% If there are cells that are still unidentifiable, 
% look at combinations of stimulation

unidentified_cells = find( unique_stim_locations==0);
% Expand the list of unidentified cells to include their neighbours 


if ~isempty(unidentified_cells) % directly stim on nulcei for other cells  

% Construct denses grids around these cells 

grid_jitters = zeros(num_per_grid_dense,2);
for i_grid = 1:num_per_grid_dense
   grid_jitters(i_grid,:)=[sin(2*pi*i_grid/num_per_grid_dense) cos(2*pi*i_grid/num_per_grid_dense)]; 
end
grid_jitters=[grid_jitters zeros(num_per_grid_dense,1)];
target_locations_dense = zeros(length(unidentified_cells)*(3*num_per_grid_dense+1),3);
for i_cell_index=1:length(unidentified_cells)
    i_cell= cell_list(unidentified_cells(i_cell_index));
    nucleus_loc=cell_locations(i_cell,:);
    grid_locs=nucleus_loc;
    grid_locs=[grid_locs;bsxfun(@plus,grid_jitters*r1,nucleus_loc)];
    grid_locs=[grid_locs;bsxfun(@plus,grid_jitters*r2,nucleus_loc)];
    grid_locs=[grid_locs;bsxfun(@plus,grid_jitters*r3,nucleus_loc)];
    target_idx=(i_cell_index-1)*(3*num_per_grid_dense+1) +(1: (3*num_per_grid_dense+1));
    target_locations_dense(target_idx,:) = grid_locs;
end
target_locations_dense(:,3)= mean(cell_locations(cell_list,3));

cell_params.locations =  cell_locations(cell_list,:);
cell_params.shape_gain = ones(length(cell_list),1);
cell_template = struct();
cell_template.shape= shape_template;
% [pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
%     cell_template,target_locations);
[pi_target_dense, ~] = get_weights_v2(cell_params, ...
    cell_template,target_locations_dense);

 stims_by_locations=zeros(size(pi_target_dense,2),size(pi_target_dense,1));
  
prob_cells_by_locations=zeros(size(pi_target_dense,2),size(pi_target_dense,1));
for l = 1:size(pi_target_dense,2)
    for i_cell_index = 1:length(cell_list)
             stims_by_locations(l,i_cell_index)=pi_target_dense(i_cell_index,l)*power_level(1);
        prob_cells_by_locations(l,i_cell_index)= stim_to_prob(...
            pi_target_dense(i_cell_index,l)*power_level(1),stim_unique,prob_trace);
    end
end

% define the neighbours of the cells ...
cell_neighbours=diag(ones(length(cell_list),1));

% A heuristic approach
% For each unidentifiable cell, find rows that maximize its firing
% probabilities subtracting the second largest one
prob_temp=prob_cells_by_locations(:,unidentified_cells);
for i_cell_index = 1:length(unidentified_cells)
    i_cell = unidentified_cells(i_cell_index);
    temp=prob_temp;
    temp(:,i_cell_index)=0;
    temp2=prob_temp(:,i_cell_index);
    temp2(temp2<1e-2)=0;
    max_diff= temp2- max(temp')';
   % pick one randomly among ties 
    [~, tar_idx]=max(max_diff);
    unique_stim_locations(i_cell)=tar_idx+size(pi_target,2);
    unique_stim_powers(i_cell)=power_level(1);
%     cell_neighbours(i_cell,prob_temp(tar_idx,:)>1e-2)=1;
    
    cell_neighbours(i_cell,stims_by_locations(tar_idx,:)>stim_threshold)=1;
end
 cell_neighbours= 1*((cell_neighbours+ cell_neighbours')>0);
 cell_neighbours = 1*(expm(cell_neighbours)>0);
unidentified_cells=find(sum(cell_neighbours,1)>1);
 target_locations_all= [target_locations; target_locations_dense];

else  
    
 target_locations_all= [target_locations];
    
end
 % Summarize the stim locations and power 
 target_locations_selected=target_locations_all(unique_stim_locations,:);
 power_selected=unique_stim_powers;
 
cell_params.locations =  cell_locations(cell_list,:);
cell_params.shape_gain = ones(length(cell_list),1);
cell_template = struct();
cell_template.shape= shape_template;
% [pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
%     cell_template,target_locations);
[pi_target_selected, inner_normalized_products] = get_weights_v2(cell_params, ...
    cell_template,target_locations_selected);


%----------------------------------------------------------------%
% Find spots on the planes of each neurons  

% Construct 3d grids around these cells 
% x=r sin(\theta) cos(\phi);y=r sin(\theta) sin(\phi); z=r cos(\theta);


cell_neighbours=diag(ones(length(cell_list),1));

if ~isempty(unidentified_cells) % directly stim on nulcei for other cells  
    grid_jitters = zeros(num_z_grid*num_per_grid_dense,3);
    for i_grid = 1:num_per_grid_dense
        for z_grid = 1:num_z_grid
            grid_jitters( (z_grid-1)*num_per_grid_dense+i_grid,:)=[sin(pi*z_grid/num_z_grid)*cos(2*pi*i_grid/num_per_grid_dense)...
                sin(pi*z_grid/num_z_grid)*sin(2*pi*i_grid/num_per_grid_dense) cos(pi*z_grid/num_z_grid)];
        end
    end
   
    grid_jitters=unique(round(grid_jitters,3),'rows');
    stim_locations_3d=zeros(length(unidentified_cells),1);
    stim_powers_3d= power_level(1)*ones(length(unidentified_cells),1);
    
    target_locations_3d = zeros(length(unidentified_cells)*(3*size(grid_jitters,1)+1),3);
    for i_cell_index=1:length(unidentified_cells)
        
        i_cell= cell_list(unidentified_cells(i_cell_index));
        nucleus_loc = cell_locations(i_cell,:);
        grid_locs=nucleus_loc;
        grid_locs=[grid_locs;bsxfun(@plus,grid_jitters*r1,nucleus_loc)];
        grid_locs=[grid_locs;bsxfun(@plus,grid_jitters*r2,nucleus_loc)];
        grid_locs=[grid_locs;bsxfun(@plus,grid_jitters*r3,nucleus_loc)];
        target_idx=(i_cell_index-1)*(3*size(grid_jitters,1)+1) +(1: (3*size(grid_jitters,1)+1));
        target_locations_3d(target_idx,:) = grid_locs;
    end
    cell_params.locations =  cell_locations(cell_list,:);
    cell_params.shape_gain = ones(length(cell_list),1);
    cell_template = struct();
    cell_template.shape= shape_template;
    [pi_target_3d, ~] = get_weights_v2(cell_params, ...
        cell_template,target_locations_3d);
    
    stims_by_locations=zeros(size(pi_target_3d,2),length(cell_list));
    for l = 1:size(pi_target_3d,2)
        for i_cell_index = 1:length(cell_list)
            stims_by_locations(l,i_cell_index)=pi_target_3d(i_cell_index,l)*power_level(1);
        end
    end
    
%     prob_cells_by_locations=zeros(size(pi_target_3d,2),length(cell_list));
%     for l = 1:size(pi_target_3d,2)
%         for i_cell_index = 1:length(cell_list)
%             prob_cells_by_locations(l,i_cell_index)= stim_to_prob(...
%                 pi_target_3d(i_cell_index,l)*power_level(1),stim_unique,prob_trace);
%         end
%     end
    

    stim_locations_3d_unidentified=zeros(length(cell_list),3);
    for i_cell_idx = 1:length(unidentified_cells)
        i_cell= unidentified_cells(i_cell_idx);
        temp=stims_by_locations;
%         temp=prob_cells_by_locations;
        temp(:,i_cell)=0;
        temp2=stims_by_locations(:,i_cell);
%         temp2(temp2<eff_stim_threshold)=0;
%         temp(temp<eff_stim_threshold)=0;
%         temp2=prob_cells_by_locations(:,i_cell);
%         temp2(temp2<1e-2)=0;
        %max_diff= (temp2- max(temp')')./( max(temp')'+1).*temp2;
        max_diff= (temp2- max(temp')');
        % find the location maximize the stim received while suppressing
        % the rest 
        %max_diff= temp2- max(temp')';
        %max_diff(max(temp')'> stim_threshold)=0;
        
        % pick one randomly among ties
        [~, tar_idx]=max(max_diff);
        stim_locations_3d_unidentified(i_cell,:)=target_locations_3d(tar_idx,:);
         cell_neighbours(i_cell,stims_by_locations(tar_idx,:)>stim_threshold)=1;
    end
     cell_neighbours= 1*((cell_neighbours+ cell_neighbours')>0);
 cell_neighbours = 1*(expm(cell_neighbours)>0);

end


% Use the heuristic approach to pick the stim spot, but this time using the
% max diff between stim size 

% A heuristic approach
% For each unidentifiable cell, find rows that maximize its firing
% probabilities subtracting the second largest one

 % Summarize the stim locations and power 
 target_locations_3d_selected=zeros(length(cell_list),3);
 for i_cell_idx = 1:length(cell_list)
     if sum(unidentified_cells==i_cell_idx)>0
        target_locations_3d_selected(i_cell_idx,:)=stim_locations_3d_unidentified(i_cell_idx,:);
     else
        target_locations_3d_selected(i_cell_idx,:)=cell_locations(cell_list(i_cell_idx),:);    
     end
 end
 power_3d_selected=power_level(1)*ones(length(cell_list),1);

cell_params.locations =  cell_locations(cell_list,:);
cell_params.shape_gain = ones(length(cell_list),1);
cell_template = struct();
cell_template.shape= shape_template;
% [pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
%     cell_template,target_locations);
[pi_target_3d, inner_normalized_products_3d] = get_weights_v2(cell_params, ...
    cell_template,target_locations_3d_selected);


end

% Check if the new stimulations are non-degenerate:
% eig(prob_cells_by_locations(unique_stim_locations(unidentified_cells)-size(pi_target,2),unidentified_cells))
% prob_cells_by_locations(unique_stim_locations(unidentified_cells)-size(pi_target,2),unidentified_cells)
% prob_cells_by_locations(unique_stim_locations(unidentified_cells),:)

% If it is still degenerate, we will pass these pairs to the full inference
% %%
% figure(1)
% scatter(target_locations{this_plane}(:,2),...
%     target_locations{this_plane}(:,1),...
%     'Marker','o','SizeData',3,...
%     'MarkerFaceColor','b', 'MarkerEdgeColor','b', 'MarkerFaceAlpha',0.2)
% 
% hold on;
% scatter(target_locations_dense{this_plane}(:,2),target_locations_dense{this_plane}(:,1),...
%     target_locations_dense{this_plane}(:,2),target_locations_dense{this_plane}(:,2),...
%     'Marker','o','SizeData',3,...
%    'MarkerFaceColor','g', 'MarkerEdgeColor','green','MarkerFaceAlpha',0.2)
% 
% hold on;
% scatter(cell_locations(cell_group_list{this_plane},2),...
%     cell_locations(cell_group_list{this_plane},1),...
%     'Marker','o','SizeData',10,...
%     'MarkerFaceColor','r', 'MarkerEdgeColor','r', 'MarkerFaceAlpha',1)
% 
% hold on;
% 
% 
% scatter(target_locations{this_plane}(unique_stim_locations(setdiff(1:end,unidentified_cells)),2),...
%     target_locations{this_plane}(unique_stim_locations(setdiff(1:end,unidentified_cells)),1),...
%     'Marker','o','SizeData',20,...
%     'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',1)
% hold on;
% scatter(target_locations_dense{this_plane}( unique_stim_locations(unidentified_cells)-size(pi_target,2),2),...
%     target_locations_dense{this_plane}( unique_stim_locations(unidentified_cells)-size(pi_target,2),1),...
%     'Marker','o','SizeData',20,...
%    'MarkerFaceColor','k', 'MarkerEdgeColor','k','MarkerFaceAlpha',1)
% hold off;
% 
% 
% xlabel('X (um)');
% ylabel('Y (um)');
% 
% axis off;
% % saveas(1,strcat('./Figures/Preprocess/','Selected_locations','.jpg'));
