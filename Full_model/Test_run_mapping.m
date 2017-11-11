addpath(genpath('../GitHub/mapping-inference'),genpath('../GitHub/odessa-beta-beta'));
%% Run get_experiment_setup:

disp('Test get_experiment_setup')
experiment_setup=get_experiment_setup();
experiment_setup.is_exp = 0;
experiment_setup.enable_user_breaks = 0;
%% Generate cells
disp('Test generate_neurons')
experiment_setup.sim=get_simulation_setup();
experiment_setup.neurons=generate_neurons(experiment_setup);
%% Initialize neighbourhoods 
disp('Test create_neighbourhoods')
neighbourhoods = create_neighbourhoods(experiment_setup);

%%
this_neighbourhood=neighbourhoods(8);
this_group='undefined';
group_profile=experiment_setup.groups.(this_group);

%manually give a batch ID to neighbourhood
this_neighbourhood.batch_ID=1;
 
%% Test designing new trials:

disp('Test design trials')
experiment_query_this_group=design_undefined(this_neighbourhood,group_profile);
%% Test simulation:
 disp('Test simulation')
 
 
 %% Test model fitting
 disp('Test inference')
 
 
% detect pscs and run vi
for i = 1:length(group_names)
    
    this_group = group_names{i};
    this_exp_query = experiment_query.(this_group);
    group_profile=experiment_setup.groups.(this_group);
    
    experiment_query.(this_group) = ...
        experiment_setup.groups.(this_group).psc_detect_function(this_exp_query,neighbourhood, group_profile, experiment_setup);
    
    neighbourhood = ...
        experiment_setup.groups.(this_group).inference_function(this_exp_query,neighbourhood,group_profile, experiment_setup.prior_info);
    
end

 %% Test regrouping
 
 disp('Test regrouping')
 
 
% regroup cells
for i = 1:length(group_names)
    %for j = setdiff(1:length(group_names),i)
    this_group = group_names{i};
    %         to_group = group_names{j};
    to_groups=setdiff(group_names,this_group);
    group_profile=experiment_setup.groups.(this_group);
    %   regroup_func = experiment_setup.groups.(this_group).regroup_functions.(to_group);
    neighbourhood  = experiment_setup.groups.(this_group).regroup_functions(neighbourhood,to_groups,group_profile);
%     end
end


 