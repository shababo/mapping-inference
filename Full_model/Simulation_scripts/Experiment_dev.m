%% One instance in the experiment:

% experiment_setup, neighbourhoods(neighbourhood_ID)

this_neighbourhood= neighbourhoods(neighbourhood_ID);
%% Grab new members

% Need a function for this 

%% Find out what groups there are 
unique_group_names=unique([this_neighbourhood.neurons(:).group_type_ID]);

% for i_group =1:length(unqiue_group_names)
    group_type_ID=unique_group_names(i_group);
    this_group_ID=index([experiment_setp.groups(:).group_type_ID]==group_type_ID); 
    group_profile=experiment_setup.groups(this_group_ID);

%% for each analysis group:
% design new trials
[experiment_query] = group_profile.design_function(this_neighbourhood,group_profile);
% calculate holograms
[experiment_query]=calculate_holograms(experiment_query);
% collect data:
switch experiment_setup.experiment_type
    case 'Simulation'
        [experiment_query] = draw_samples(experiment_query,this_neighbourhood);
    otherwise
end
% send it through OASIS
[experiment_query]=event_detection(experiment_query);
% record new trials:
[experiment_records]=record_experiments(experiment_records,experiment_query);

% Analysis:
[this_neighbourhood]=group_profile.inference_function(experiment_records,this_neighbourhood,group_profile);
    
 [this_neighbourhood]=group_profile.regroup_function(this_neighbourhood,group_profile);
    
