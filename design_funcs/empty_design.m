function experiment_query = empty_design(neighbourhood, group_profile)

group_ID=group_profile.group_ID;
cells_this_group=index([this_neighbourhood.neurons(:).group_ID]==group_ID);
number_cells_this_group=length(cells_this_group);
number_cells_all= length(this_neighbourhood.neurons);
loc_counts=zeros(number_cells_this_group,1);

experiment_query.neighbourhood_ID = neighbourhood.neighbourhood_ID;
experiment_query.data_gen_type = 'acq';
experiment_query.trials = [];
experiment_query.batch_info.batch_id = 1;
experiment_query.batch_info.trial_counts = 0;
experiment_query.batch_info.duration = NaN;
experiment_query.batch_info.hologram_computing_time = NaN;