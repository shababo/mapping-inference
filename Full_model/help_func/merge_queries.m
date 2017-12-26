function [full_queries]=merge_queries(neigbhourhood_ID,experiment_queries)

full_queries=struct;
full_queries.undefined=struct;
full_queries.connected=struct;
for j = 1:length(neigbhourhood_ID)
    i_neighbourhood = neigbhourhood_ID(j);
    for i = 2:size(experiment_queries,2)
        if (isfield(experiment_queries(i_neighbourhood,i).undefined,'trials'))
            if (~isfield(full_queries.undefined,'trials'))
                full_queries.undefined.trials=experiment_queries(i_neighbourhood,i).undefined.trials;
            else
                n_exist = length(full_queries.undefined.trials);
                n_new = length(experiment_queries(i_neighbourhood,i).undefined.trials);
                full_queries.undefined.trials(n_exist+(1:n_new))=experiment_queries(i_neighbourhood,i).undefined.trials;
            end
        end
        
        if (isfield(experiment_queries(i_neighbourhood,i).connected,'trials'))
            if (~isfield(full_queries.connected,'trials'))
                full_queries.connected.trials=experiment_queries(i_neighbourhood,i).connected.trials;
            else
                n_exist = length(full_queries.connected.trials);
                n_new = length(experiment_queries(i_neighbourhood,i).connected.trials);
                full_queries.connected.trials(n_exist+(1:n_new))=experiment_queries(i_neighbourhood,i).connected.trials;
            end
        end
    end
end



