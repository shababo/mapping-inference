function [experiment_queries, neighbourhoods] = collect_exp_queries(experiment_setup,neighbourhood_IDs,batch_IDs,varargin)

if ~isempty(varargin) && ~isempty(varargin{1})
    
    dirpath = varargin{1};
 
else
    
    dirpath = experiment_setup.analysis_root;
end

group_names = experiment_setup.group_names;

for j = 1:length(neighbourhood_IDs)
    
    for i = 1:length(batch_IDs)

        batchsavepath = [dirpath experiment_setup.exp_id ...
                        '_n' num2str(neighbourhood_IDs(j))...
                        '_b' num2str(batch_IDs(i)) '_complete.mat'];

        load(batchsavepath)
        
        if i == 1
            for k = 1:length(group_names)
                if ~isfield(experiment_query,group_names{k})
                    experiment_query.(group_names{k}) = struct();
                end
            end
            experiment_query.batch_trial_rate = 0;
        end
        
        neighbourhoods(j,i) = neighbourhood;
        experiment_queries(j,i) = experiment_query;


    end
end


