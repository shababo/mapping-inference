function neighbourhoods = create_neighbourhoods(experiment_setup)
% initialize the neighbourhoods 
 
cell_locations=reshape([experiment_setup.neurons.location],length(experiment_setup.neurons(1).location),[])';
% assignin('base','cell_locations',cell_locations)
% assignin('base','experiment_setup',experiment_setup)
z_min=min(cell_locations(:,3));
z_max=max(cell_locations(:,3));
z_slide_width=experiment_setup.neighbourhood_params.height;
z_borders=[z_min:z_slide_width:z_max];% z_max];
% assignin('base','z_borders',z_borders)
number_of_neighbourhoods = length(z_borders)-1;

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
        
   location_vec = [neighbourhoods(i_neighbourhood).neurons.location];
   z_locs = location_vec(3:3:end);
   neighbourhoods(i_neighbourhood).center = [0 0 mean(z_locs)]; % mean(z_borders(i_neighbourhood:i_neighbourhood+1))];
   neighbourhoods(i_neighbourhood).computing_time=struct;
   neighbourhoods(i_neighbourhood).batch_ID=1;
   
end


% % Include nearby cells 
% 
for i_neighbourhood = 1:number_of_neighbourhoods
   %neighbourhoods(i_neighbourhood)=struct;
   nearby_cell_list = find(cell_locations(:,3)<  neighbourhoods(i_neighbourhood).center(3)+experiment_setup.neighbourhood_params.buffer_height/2 & ...
       cell_locations(:,3)>  neighbourhoods(i_neighbourhood).center(3)-experiment_setup.neighbourhood_params.buffer_height/2);
   secondary_cell_list= setdiff(nearby_cell_list, [neighbourhoods(i_neighbourhood).neurons(:).cell_ID]);
   number_of_prim_cells =length(neighbourhoods(i_neighbourhood).neurons);
   for i_cell = 1:length(secondary_cell_list)
       cell_ID=secondary_cell_list(i_cell);
       neighbourhoods(i_neighbourhood).neurons(i_cell+number_of_prim_cells)=experiment_setup.neurons(cell_ID);
   end
   for i_cell = 1:length(neighbourhoods(i_neighbourhood).neurons)
       if i_cell > number_of_prim_cells
           neighbourhoods(i_neighbourhood).neurons(i_cell).group_ID='secondary';
       else
           neighbourhoods(i_neighbourhood).neurons(i_cell).group_ID=experiment_setup.default_group; % all are initialized as undefined
       end
   end
end



% Calculate the candidate grid for each cell in each neigbourhood
group_names = experiment_setup.group_names;
for i_neighbourhood = 1:number_of_neighbourhoods
    for i_cell = 1:length(neighbourhoods(i_neighbourhood).neurons)
        neighbourhoods(i_neighbourhood).neurons(i_cell).stim_locations=struct([]);
        neighbourhoods(i_neighbourhood).neurons(i_cell).group_ID=experiment_setup.default_group; % all are initialized as undefined
        neighbourhoods(i_neighbourhood).neurons(i_cell).primary_indicator=true;
        neighbourhoods(i_neighbourhood).neurons(i_cell).stim_locations=struct;
        for i_group = 1:length(group_names) % if this is costly we could move to when a neuron is initialized into a group
            if isfield(experiment_setup.groups.(group_names{i_group}),'design_func_params')
                design_params=experiment_setup.groups.(group_names{i_group}).design_func_params;
                cell_params=neighbourhoods(i_neighbourhood).neurons;
                neighbourhoods(i_neighbourhood).neurons(i_cell).stim_locations.(group_names{i_group})=...
                    get_stim_locations(i_cell,group_names{i_group},cell_params,design_params,experiment_setup.prior_info.template_cell,...
                                        neighbourhoods(i_neighbourhood).center(3),experiment_setup.exp.foe_bounds);
%                 if experiment_setup.is_exp
%                     these_locs = neighbourhoods(i_neighbourhood).neurons(i_cell).stim_locations.(group_names{i_group}).grid;
%                     for i_dim = 1:size(experiment_setup.exp.foe_bounds,1)
%                         these_locs(these_locs(:,i_dim) < experiment_setup.exp.foe_bounds(i_dim,1),i_dim) = experiment_setup.exp.foe_bounds(i_dim,1);
%                         these_locs(these_locs(:,i_dim) > experiment_setup.exp.foe_bounds(i_dim,2),i_dim) = experiment_setup.exp.foe_bounds(i_dim,2);
%                     end
%                     neighbourhoods(i_neighbourhood).neurons(i_cell).stim_locations.(group_names{i_group}).grid = these_locs;    
%                 end
            end
        end
    end
end


% Initialize PR and gain parameters 
for i_neighbourhood = 1:number_of_neighbourhoods
    for i_cell = 1:length(neighbourhoods(i_neighbourhood).neurons)
        current_params=experiment_setup.prior_info.PR_prior;
        group_profile=experiment_setup.groups.(neighbourhoods(i_neighbourhood).neurons(i_cell).group_ID);
        bounds= group_profile.inference_params.bounds.PR;
        quantile_prob=group_profile.regroup_func_params.quantile_prob;
        neighbourhoods(i_neighbourhood).neurons(i_cell).PR_params(1)=calculate_posterior(...
            current_params,bounds,quantile_prob);

        current_params=experiment_setup.prior_info.gain_prior;
        group_profile=experiment_setup.groups.(neighbourhoods(i_neighbourhood).neurons(i_cell).group_ID);
        bounds= group_profile.inference_params.bounds.gain;
        quantile_prob=group_profile.regroup_func_params.quantile_prob;
        neighbourhoods(i_neighbourhood).neurons(i_cell).gain_params(1)=calculate_posterior(...
            current_params,bounds,quantile_prob);

    end
end
%                 neighbourhoods(i_neighbourhood).neurons(i_cell).PR_params(1).mean=0;
%                 neighbourhoods(i_neighbourhood).neurons(i_cell).PR_params(1).variance=0;
%                 neighbourhoods(i_neighbourhood).neurons(i_cell).PR_params(1).upper_quantile=0;
%                 neighbourhoods(i_neighbourhood).neurons(i_cell).PR_params(1).lower_quantile=0;
%                 neighbourhoods(i_neighbourhood).neurons(i_cell).PR_params(1).nonzero_prob=0;
                
                %neighbourhoods(i_neighbourhood).neurons(i_cell).gain=struct([]);
         

