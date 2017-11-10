function [experiment_query] = generate_psc_data(experiment_query,experiment_setup,this_neighbourhood)

group_names= fieldnames(experiment_query);
for i_group = 1:length(group_names)
    experiment_query_this_group =experiment_query.(group_names(i_group));
    experiment_query_this_group = draw_point_processes(experiment_query_this_group, ...
        this_neighbourhood,experiment_setup);
    experiment_query.(group_names(i_group))=experiment_query_this_group;
end

end



