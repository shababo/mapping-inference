function [neurons,stimuli_size] = get_stim_size_batch(neurons,trials,simulation_params)
%% Calculate stimulation size for simulation
min_dist = 1;
number_of_trials = length(trials);
number_of_cells =length(neurons);

stimuli_size=zeros(number_of_trials,number_of_cells);
for i_cell = 1:number_of_cells
    neurons(i_cell).truth.shapes_sim=struct;
    neurons(i_cell).truth.shapes_sim.locations=zeros(0,3);
    neurons(i_cell).truth.shapes_sim.values=[];
    neurons(i_cell).truth.shapes_sim.indices=cell([0 1]);
end

for l=1:number_of_trials
    this_loc=trials(l).locations; % 1 spot per trial
    this_power=trials(l).power_levels; % 1 spot per trial
    for i_cell=1:number_of_cells
        rel_pos =  this_loc - neurons(i_cell).truth.location;
        
        existing_loc_dim=size(neurons(i_cell).truth.shapes_sim.locations);
        sq_dist=sum((  neurons(i_cell).truth.shapes_sim.locations- ones(existing_loc_dim(1),1)*rel_pos).^2,2);
        if (existing_loc_dim(1) == 0)   | min(sq_dist.^(1/2))> min_dist
            tmp = [];
        else
            [~,tmp]=min(sq_dist.^(1/2));
        end
        
        if isempty(tmp) % this is a new location:
            neurons(i_cell).truth.shapes_sim.locations= [neurons(i_cell).truth.shapes_sim.locations; rel_pos];
            this_size = griddata(simulation_params.mesh_grid(:,1),simulation_params.mesh_grid(:,2),simulation_params.mesh_grid(:,3),...
                neurons(i_cell).truth.shape,rel_pos(1),rel_pos(2),rel_pos(3),'linear');
            this_size = neurons(i_cell).truth.excite(this_size);
            neurons(i_cell).truth.shapes_sim.values=[ neurons(i_cell).truth.shapes_sim.values; this_size];
            neurons(i_cell).truth.shapes_sim.indices{length(neurons(i_cell).truth.shapes_sim.indices)+1}=l;
        else
            this_size = neurons(i_cell).truth.shapes_sim.values(tmp);
            neurons(i_cell).truth.shapes_sim.indices{tmp}=[neurons(i_cell).truth.shapes_sim.indices{tmp}; l];
        end
        
        stimuli_size(l,i_cell) = this_power*this_size;
    end
    
    if mod(l,50)==1
        fprintf('Trials evaluated: %d;\n',l)
    end
end


%------------------------------------------------------------------------%
% Use get_stim_size_batch_prior to predict the stim size given prior distribution
% Use get_stim_size_batch_expected to predict the stim size given current posterior
