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
experiment_query_this_group=design_undefined(this_neighbourhood,group_profile,experiment_setup);

%% 
experiment_query.undefined=experiment_query_this_group;
%%
experiment_setup.patched_neuron=struct;
    experiment_setup.patched_neuron.background_rate=1e-4;
    experiment_setup.patched_neuron.cell_type=[];
%% Test simulation:
 disp('Test simulation')
 
 [experiment_query] = generate_psc_data(experiment_query,experiment_setup,this_neighbourhood);

 %% manually move event times to this_experiment_query
 n_trials = length(experiment_query.undefined.trials);
 
 for i_trial = 1:n_trials
    experiment_query.undefined.trials(i_trial).event_times=experiment_query.undefined.trials(i_trial).truth.event_times;
    
 end
 
 %%
 experiment_query_this_group =experiment_query.undefined;
 %% Test model fitting
 disp('Test inference')
 prior_info=experiment_setup.prior_info;
  this_neighbourhood=inference_undefined(experiment_query_this_group,this_neighbourhood,group_profile,experiment_setup);
 
 %% Test regrouping
 
 disp('Test regrouping')
 
 this_group = 'undefined';
 i=1;
 group_names={'undefined','disconnected','connected','alive'};   
 group_profile=get_undefined();
 %%
 for j = setdiff(1:length(group_names),i)
     to_group=group_names{j};
     if isfield( group_profile.regroup_functions, to_group)
     regroup_func = group_profile.regroup_functions.(to_group);
     this_neighbourhood  = regroup_func(this_neighbourhood,group_profile);
     %     end
     end
 end


 