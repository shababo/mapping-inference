function experiment_query = empty_design(neighbourhood)

experiment_query.neighbourhood_id = neighbourhood.neighbourhood_id;
experiment_query.data_gen_type = 'acq';
experiment_query.trials = [];
experiment_query.batch_info.batch_id = 1;
experiment_query.batch_info.trial_counts = 0;
experiment_query.batch_info.duration = NaN;
experiment_query.batch_info.hologram_computing_time = NaN;