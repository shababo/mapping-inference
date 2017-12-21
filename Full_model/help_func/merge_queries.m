function [full_queries]=merge_queries(neigbhourhood_ID,experiment_queries)

full_queries=struct;
full_queries.undefined=struct;
full_queries.connected=struct;

for i = 2:size(experiment_queries,2)
    if (isfield(experiment_queries(neigbhourhood_ID,i).undefined,'trials'))
        if (~isfield(full_queries.undefined,'trials'))
            full_queries.undefined.trials=experiment_queries(neigbhourhood_ID,i).undefined.trials;
        else
            n_exist = length(full_queries.undefined.trials);
            n_new = length(experiment_queries(neigbhourhood_ID,i).undefined.trials);
            full_queries.undefined.trials(n_exist+(1:n_new))=experiment_queries(neigbhourhood_ID,i).undefined.trials; 
        end
    end
    
    if (isfield(experiment_queries(neigbhourhood_ID,i).connected,'trials'))
        if (~isfield(full_queries.connected,'trials'))
            full_queries.connected.trials=experiment_queries(neigbhourhood_ID,i).connected.trials;
        else
            n_exist = length(full_queries.connected.trials);
            n_new = length(experiment_queries(neigbhourhood_ID,i).connected.trials);
            full_queries.connected.trials(n_exist+(1:n_new))=experiment_queries(neigbhourhood_ID,i).connected.trials; 
        end
    end
end 



