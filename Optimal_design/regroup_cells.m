function data = regroup_cells(data)

i = data.design.i;
params = data.params;

undefined_to_disconnected = ...
    intersect(find(data.design.mean_gamma_undefined<params.design.disconnected_threshold),find( data.design.undefined_cells{i}{data.design.iter}));
undefined_to_connected = ...
    intersect(find(data.design.mean_gamma_undefined>params.design.connected_threshold),find( data.design.undefined_cells{i}{data.design.iter}));
% cells move together with their neighbours
%             undefined_to_disconnected=find(sum(cell_neighbours(undefined_to_disconnected,:),1)>0)';
%             undefined_to_connected =find(sum(cell_neighbours(undefined_to_connected,:),1)>0);
%             % if there are conflicts, move them to the potentially connected cells
%             undefined_to_disconnected=setdiff(undefined_to_disconnected,undefined_to_connected);

disconnected_to_undefined = intersect(find(data.design.mean_gamma_disconnected>params.design.disconnected_confirm_threshold),...
    find(data.design.potentially_disconnected_cells{i}{data.design.iter}));
disconnected_to_dead = intersect(find(data.design.mean_gamma_disconnected<params.design.disconnected_confirm_threshold),...
    find(data.design.potentially_disconnected_cells{i}{data.design.iter}));

%             disconnected_to_undefined=find(sum(cell_neighbours(disconnected_to_undefined,:),1)>0);
%             % if there are conflicts, move them to the potentially connected cells
%             disconnected_to_dead=setdiff(disconnected_to_dead,disconnected_to_undefined);


connected_to_dead = intersect(find(data.design.mean_gamma_connected<params.design.disconnected_confirm_threshold),...
    find(data.design.potentially_connected_cells{i}{data.design.iter}));
connected_to_alive = intersect(find(data.design.mean_gamma_connected>params.design.connected_confirm_threshold),...
    find(data.design.potentially_connected_cells{i}{data.design.iter}));
%             change_gamma =abs(data.design.gamma_path{i}(:,data.design.iter+1)-data.design.gamma_path{i}(:,data.design.iter));
connected_to_alive = intersect(find(data.design.change_gamma<params.design.change_threshold),...
    connected_to_alive);

% Eliminate the weakly identifiable pairs if they are both assign to a
% group:
%moved_cells = [connected_to_dead; connected_to_alive]';
%cells_and_neighbours=find(sum(cell_neighbours(moved_cells,:),1)>0);
%neighbours_not_included=intersect(find(data.design.potentially_connected_cells{i}{data.design.iter}), setdiff(cells_and_neighbours,moved_cells));
%blacklist=find(sum(cell_neighbours(neighbours_not_included,:),1)>0);
%connected_to_dead=setdiff(connected_to_dead ,blacklist);
%connected_to_alive=setdiff(connected_to_alive,blacklist);

% Update the cell lists:
data.design.undefined_cells{i}{data.design.iter+1}=data.design.undefined_cells{i}{data.design.iter};
data.design.undefined_cells{i}{data.design.iter+1}(undefined_to_disconnected)=0;data.design.undefined_cells{i}{data.design.iter+1}(undefined_to_connected)=0;
data.design.undefined_cells{i}{data.design.iter+1}(disconnected_to_undefined)=1;

data.design.potentially_disconnected_cells{i}{data.design.iter+1}=data.design.potentially_disconnected_cells{i}{data.design.iter};
data.design.potentially_disconnected_cells{i}{data.design.iter+1}(disconnected_to_dead)=0;data.design.potentially_disconnected_cells{i}{data.design.iter+1}(disconnected_to_undefined)=0;
data.design.potentially_disconnected_cells{i}{data.design.iter+1}(undefined_to_disconnected)=1;


data.design.potentially_connected_cells{i}{data.design.iter+1}=data.design.potentially_connected_cells{i}{data.design.iter};
data.design.potentially_connected_cells{i}{data.design.iter+1}(connected_to_dead)=0;data.design.potentially_connected_cells{i}{data.design.iter+1}(connected_to_alive)=0;
data.design.potentially_connected_cells{i}{data.design.iter+1}(undefined_to_connected)=1;

data.design.dead_cells{i}{data.design.iter+1}=data.design.dead_cells{i}{data.design.iter};
data.design.dead_cells{i}{data.design.iter+1}(disconnected_to_dead)=1;data.design.dead_cells{i}{data.design.iter+1}(connected_to_dead)=1;

data.design.alive_cells{i}{data.design.iter+1}=data.design.alive_cells{i}{data.design.iter};
data.design.alive_cells{i}{data.design.iter+1}(connected_to_alive)=1;