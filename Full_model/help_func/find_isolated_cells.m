function [cell_IDs] = find_isolated_cells(experiment_setup,isolation_zone)
%%

cell_locations=reshape([experiment_setup.neurons(:).location],3,[])';


cell_IDs=[];
flags=zeros(length(experiment_setup.neurons),3);
for i_cell = 1:length(experiment_setup.neurons)
    
    this_location = experiment_setup.neurons(i_cell).location;
    relative_dist=abs(cell_locations-ones(length(experiment_setup.neurons),1)* this_location);
    for i_dim = 1:length(isolation_zone)
        flags(:,i_dim)=relative_dist(:,i_dim)<isolation_zone(i_dim);
    end
    flags(i_cell,:)=0;
    if max(prod(flags,2))==0
       cell_IDs=[cell_IDs  experiment_setup.neurons(i_cell).cell_ID];
    end
end