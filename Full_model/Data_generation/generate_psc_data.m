function [experiment_query] = generate_psc_data(experiment_query,...
    experiment_setup,neighbourhood)
%%
group_names= fieldnames(experiment_query);
for i_group = 1:length(group_names)
    
    if any(get_group_inds(neighbourhood,group_names{i_group})) &&   isfield(experiment_query.(group_names{i_group}),'trials')
        
        if ~isempty(experiment_query.(group_names{i_group}).trials)
            experiment_query_this_group =experiment_query.(group_names{i_group});
            
            [experiment_query_this_group] = draw_point_processes_intensity(experiment_query_this_group, ...
                neighbourhood,experiment_setup);
            if experiment_setup.sim.sim_vclamp
                experiment_query_this_group = gen_trace_from_mpp(experiment_query_this_group, ...
                    neighbourhood,experiment_setup);
            end
            experiment_query.(group_names{i_group})=experiment_query_this_group;
        end
    end
end

end



